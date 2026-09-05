#!/usr/bin/env julia

using NCDatasets
using Printf
using SBLToolkit

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_INPUT = joinpath(REPO_ROOT, ".data", "gabs3", "gabls3_scm_cabauw_obs_v33.nc")
const RHO_CP_AIR = 1.225 * 1004.67

function required_variable(ds::NCDataset, name::String)
    haskey(ds, name) || throw(ArgumentError("required variable '$name' is absent from $(ds.path)"))
    return ds[name][:]
end

function to_kelvin(temperature::AbstractArray{T}, units::AbstractString) where {T<:AbstractFloat}
    lowercase(units) in ("celsius", "degc", "degree_celsius") && return temperature .+ T(273.15)
    return temperature
end

function profile_matrix(values, n_levels::Int, n_times::Int, name::String)
    length(values) == n_levels * n_times || throw(ArgumentError(
        "$name has $(length(values)) values; expected $(n_levels * n_times)",
    ))
    return reshape(Float64.(values), n_levels, n_times)
end

function load_cabauw_prescribed_flux(path::String)
    isfile(path) || throw(ArgumentError(
        "Cabauw forcing file not found at '$path'. This Level 4.5 dry-run does not synthesize observations.",
    ))
    return NCDataset(path, "r") do ds
        z = Float64.(required_variable(ds, "z"))
        time = Float64.(required_variable(ds, "time"))
        u = profile_matrix(required_variable(ds, "u"), length(z), length(time), "u")
        v = profile_matrix(required_variable(ds, "v"), length(z), length(time), "v")
        temperature_variable = ds["temp"]
        theta = to_kelvin(profile_matrix(required_variable(ds, "temp"), length(z), length(time), "temp"),
            string(get(temperature_variable.attrib, "units", "")))
        sensible_heat = profile_matrix(required_variable(ds, "hs"), length(z), length(time), "hs")

        expected_shape = (length(z), length(time))
        all(size(profile) == expected_shape for profile in (u, v, theta, sensible_heat)) ||
            throw(ArgumentError("u, v, temp, and hs must have shape $expected_shape"))
        all(isfinite, (z..., time..., u..., v..., theta..., sensible_heat...)) ||
            throw(ArgumentError("Cabauw forcing contains nonfinite values"))
        issorted(z) && all(diff(z) .> 0.0) || throw(ArgumentError("Cabauw heights must be strictly increasing"))
        issorted(time) && all(diff(time) .> 0.0) || throw(ArgumentError("Cabauw times must be strictly increasing"))
        return (; z, time, u, v, theta, sensible_heat)
    end
end

function run_prescribed_flux_dry_run(path::String=DEFAULT_INPUT)
    forcing = load_cabauw_prescribed_flux(path)
    state = initialize_cabauw_state(forcing.z)
    state.u .= forcing.u[:, 1]
    state.v .= forcing.v[:, 1]
    state.theta_v .= forcing.theta[:, 1]
    state.skin_temperature = forcing.theta[1, 1]
    gating_params = BifurcationGatingParams(
        epsilon_on=-0.10, epsilon_off=-0.05, zeta_z_tol=1e-2, km_floor=0.1, kh_floor=0.015,
    )
    gs_config = GSPTModelConfig()
    diagnostics = SCMDiagnostics{Float64}[]

    for index in 1:(length(forcing.time) - 1)
        dt = (forcing.time[index + 1] - forcing.time[index]) * 86_400.0
        theta_flux = forcing.sensible_heat[1, index] / RHO_CP_AIR
        config = SCMConfig(forcing.z; dt, theta_flux, surface_mode=:prescribed)

        # This fixture lacks TKE, radiation, and geostrophic winds. Each observed
        # profile is therefore diagnostic input; missing fields are not inferred.
        state.u .= forcing.u[:, index]
        state.v .= forcing.v[:, index]
        state.theta_v .= forcing.theta[:, index]
        push!(diagnostics, step_scm!(state, config, gating_params, gs_config; enable_gating=true))
    end
    return forcing, diagnostics
end

function print_summary(path::String, forcing, diagnostics)
    total_gated_levels = sum(diagnostic.gated_levels for diagnostic in diagnostics)
    @printf("Level 4.5 prescribed-flux ingestion dry-run\n")
    @printf("  input: %s\n", path)
    @printf("  records advanced: %d\n", length(diagnostics))
    @printf("  height levels: %d\n", length(forcing.z))
    @printf("  input surface H: %.3f to %.3f W m^-2\n",
        minimum(forcing.sensible_heat[1, :]), maximum(forcing.sensible_heat[1, :]))
    @printf("  final diagnosed T_2m: %.3f K\n", diagnostics[end].t2m)
    @printf("  gated level-timesteps: %d\n", total_gated_levels)
    println("  status: ingestion and prescribed-flux interface verified; not an empirical cold-bias benchmark")
end

input_path = isempty(ARGS) ? DEFAULT_INPUT : ARGS[1]
forcing, diagnostics = run_prescribed_flux_dry_run(input_path)
print_summary(input_path, forcing, diagnostics)