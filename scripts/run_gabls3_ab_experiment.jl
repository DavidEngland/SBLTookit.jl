using Dates
using NCDatasets
using Printf
using Statistics
using SBLToolkit

const HOURS = 24
const SECONDS_PER_HOUR = 3600
const MODEL_DT = 60.0
const DEFAULT_CABAUW_REFERENCE = joinpath(@__DIR__, "..", "data", "raw", "gabls3",
    "gabls3_scm_cabauw_obs_v33.nc")

function run_experiment(closure_mode::Symbol)
    z = collect(5.0:5.0:200.0)
    state = initialize_cabauw_state(z)
    config = SCMConfig(z; dt=MODEL_DT, surface_mode=:slab, radiation_max=400.0,
        longwave_loss=80.0, radiation_peak_time=43_200.0, soil_temperature=279.0,
        ground_conductance=2.0)
    gating_params = BifurcationGatingParams(
        epsilon_on=-0.10, epsilon_off=-0.05, zeta_z_tol=1e-2, km_floor=0.1, kh_floor=0.01,
    )
    gs_config = GSPTModelConfig()
    diagnostics = SCMDiagnostics{Float64}[]

    for _ in 1:(HOURS * SECONDS_PER_HOUR ÷ round(Int, config.dt))
        push!(diagnostics, step_scm!(state, config, gating_params, gs_config;
            closure_mode, enable_gating=closure_mode !== :control))
    end
    return diagnostics
end

function load_cabauw_t2m_reference(path::String)
    isfile(path) || throw(ArgumentError("Cabauw reference file not found at '$path'"))
    return NCDataset(path, "r") do ds
        haskey(ds, "time") || throw(ArgumentError("Cabauw reference lacks time(time)"))
        haskey(ds, "t2m") || throw(ArgumentError("Cabauw reference lacks t2m(time)"))
        time_values = ds["time"][:]
        time_hours = if time_values isa AbstractVector{<:Dates.AbstractDateTime}
            Float64.(Dates.value.(time_values .- first(time_values))) ./ 3_600_000.0
        else
            Float64.(time_values)
        end
        t2m = Float64.(ds["t2m"][:])
        length(time_hours) == length(t2m) || throw(ArgumentError("time and t2m lengths differ"))
        length(t2m) >= 2 || throw(ArgumentError("Cabauw reference requires at least two T_2m samples"))
        all(isfinite, time_hours) && all(isfinite, t2m) || throw(ArgumentError(
            "Cabauw T_2m reference contains nonfinite values",
        ))
        issorted(time_hours) && all(diff(time_hours) .> 0.0) || throw(ArgumentError(
            "Cabauw T_2m timestamps must be strictly increasing",
        ))
        units = string(get(ds["t2m"].attrib, "units", ""))
        units == "K" || throw(ArgumentError("Cabauw t2m must be in K; found '$units'"))
        return (; time_hours, t2m)
    end
end

function aligned_t2m(diagnostics::Vector{SCMDiagnostics{Float64}}, time_hours::Vector{Float64})
    all(0.0 .<= time_hours .<= HOURS) || throw(ArgumentError(
        "Cabauw timestamps lie outside the synthetic 24-hour model window",
    ))
    indices = clamp.(round.(Int, time_hours .* SECONDS_PER_HOUR ./ MODEL_DT),
        1, length(diagnostics))
    return getfield.(diagnostics[indices], :t2m)
end

function summarize_reference(label::String, diagnostics::Vector{SCMDiagnostics{Float64}}, reference)
    model_t2m = aligned_t2m(diagnostics, reference.time_hours)
    error = model_t2m .- reference.t2m
    @printf("  %s mean T_2m bias: %.3f K\n", label, mean(error))
    @printf("  %s T_2m RMSE: %.3f K\n", label, sqrt(mean(error .^ 2)))
    @printf("  %s nocturnal-minimum bias: %.3f K\n", label,
        minimum(model_t2m) - minimum(reference.t2m))
end

function summarize(label::String, diagnostics::Vector{SCMDiagnostics{Float64}})
    t2m = getfield.(diagnostics, :t2m)
    z_i = getfield.(diagnostics, :boundary_layer_height)
    surface_flux = getfield.(diagnostics, :surface_sensible_heat_flux)
    radiation = getfield.(diagnostics, :net_radiation)
    ground_flux = getfield.(diagnostics, :ground_heat_flux)
    gated_steps = sum(getfield.(diagnostics, :gated_levels))
    override_steps = sum(count.(getfield.(diagnostics, :override_mask)))
    @printf("%s\n", label)
    @printf("  final T_2m: %.3f K\n", t2m[end])
    @printf("  mean T_2m: %.3f K\n", mean(t2m))
    @printf("  mean z_i: %.2f m\n", mean(z_i))
    @printf("  mean H_0: %.2f W m^-2\n", mean(surface_flux))
    @printf("  mean R_n: %.2f W m^-2\n", mean(radiation))
    @printf("  mean G_s: %.2f W m^-2\n", mean(ground_flux))
    @printf("  gated level-timesteps: %d\n", gated_steps)
    @printf("  eligible Ri-shutdown overrides: %d\n", override_steps)
    return t2m[end]
end

control = run_experiment(:control)
gate_only = run_experiment(:gate_only)
full_gspt = run_experiment(:full_gspt)
control_t2m = summarize("Experiment A: control", control)
gate_only_t2m = summarize("Experiment B: gate-only", gate_only)
full_gspt_t2m = summarize("Experiment C: full GSPT", full_gspt)
@printf("Final T_2m difference (gate-only - control): %.3f K\n", gate_only_t2m - control_t2m)
@printf("Final T_2m difference (full GSPT - control): %.3f K\n", full_gspt_t2m - control_t2m)

if "--cabauw-reference" in ARGS
    reference_path = length(ARGS) > 1 ? ARGS[2] : DEFAULT_CABAUW_REFERENCE
    reference = load_cabauw_t2m_reference(reference_path)
    println("Observed Cabauw T_2m reference comparison")
    summarize_reference("control", control, reference)
    summarize_reference("gate-only", gate_only, reference)
    summarize_reference("full GSPT", full_gspt, reference)
    println("Reference-only calibration: the SCM remains analytically forced; this is not yet a fully observed-forcing replay.")
end