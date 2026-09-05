module GABLS3Stepper

using LinearAlgebra
using ..SBLGating
using ..GABLS3Adapters

export SCMState, SCMConfig, SCMDiagnostics, net_radiation,
    surface_sensible_heat_flux, solve_tridiagonal!, step_scm!, initialize_cabauw_state

const GRAVITY = 9.81
const RHO_CP_AIR = 1.225 * 1004.67

mutable struct SCMState{T<:AbstractFloat}
    z::Vector{T}
    u::Vector{T}
    v::Vector{T}
    theta_v::Vector{T}
    tke::Vector{T}
    skin_temperature::T
    elapsed_seconds::T
    gating_states::Vector{GatingState{T}}
end

function SCMState(
    z::AbstractVector{T},
    u::AbstractVector{T},
    v::AbstractVector{T},
    theta_v::AbstractVector{T},
    tke::AbstractVector{T},
    skin_temperature::T,
) where {T<:AbstractFloat}
    n_levels = length(z)
    n_levels >= 2 || throw(ArgumentError("SCM state requires at least two vertical levels"))
    all(length(profile) == n_levels for profile in (u, v, theta_v, tke)) || throw(ArgumentError(
        "SCM state profiles must share the vertical-grid length",
    ))
    all(isfinite, (z..., u..., v..., theta_v..., tke..., skin_temperature)) || throw(ArgumentError(
        "SCM state values must be finite",
    ))
    issorted(z) && all(diff(z) .> zero(T)) || throw(ArgumentError(
        "SCM vertical levels must be strictly increasing",
    ))
    return SCMState(
        collect(z), collect(u), collect(v), collect(theta_v), max.(collect(tke), zero(T)),
        skin_temperature, zero(T), [GatingState(T) for _ in 1:n_levels],
    )
end

struct SCMConfig{T<:AbstractFloat}
    dt::T
    dz::Vector{T}
    z0::T
    theta_flux::T
    surface_radiative_tendency::T
    surface_heat_capacity::T
    km_background::T
    kh_background::T
    geostrophic_u::T
    geostrophic_v::T
    relaxation_time::T
    bulk_richardson_critical::T
    surface_mode::Symbol
    radiation_max::T
    longwave_loss::T
    radiation_peak_time::T
    soil_temperature::T
    ground_conductance::T
end

function SCMConfig(
    z::AbstractVector{T};
    dt::Real=60.0,
    z0::Real=0.15,
    theta_flux::Real=-0.02,
    surface_radiative_tendency::Real=-2e-4,
    surface_heat_capacity::Real=2e5,
    km_background::Real=0.01,
    kh_background::Real=0.001,
    geostrophic_u::Real=8.0,
    geostrophic_v::Real=0.0,
    relaxation_time::Real=21_600.0,
    bulk_richardson_critical::Real=0.25,
    surface_mode::Symbol=:prescribed,
    radiation_max::Real=400.0,
    longwave_loss::Real=80.0,
    radiation_peak_time::Real=43_200.0,
    soil_temperature::Real=279.0,
    ground_conductance::Real=2.0,
) where {T<:AbstractFloat}
    length(z) >= 2 && issorted(z) && all(diff(z) .> zero(T)) || throw(ArgumentError(
        "SCM configuration requires strictly increasing levels",
    ))
    values = T.((dt, z0, theta_flux, surface_radiative_tendency, surface_heat_capacity,
        km_background, kh_background, geostrophic_u, geostrophic_v, relaxation_time,
        bulk_richardson_critical, radiation_max, longwave_loss, radiation_peak_time,
        soil_temperature, ground_conductance))
    all(isfinite, values) || throw(ArgumentError("SCM configuration values must be finite"))
    dt_t, z0_t, theta_flux_t, radiative_tendency_t, heat_capacity_t, km_t, kh_t,
    geostrophic_u_t, geostrophic_v_t, relaxation_time_t, critical_ri_t, radiation_max_t,
    longwave_loss_t, radiation_peak_time_t, soil_temperature_t, ground_conductance_t = values
    dt_t > zero(T) || throw(ArgumentError("dt must be positive"))
    z0_t > zero(T) || throw(ArgumentError("z0 must be positive"))
    heat_capacity_t > zero(T) || throw(ArgumentError("surface_heat_capacity must be positive"))
    km_t >= zero(T) && kh_t >= zero(T) || throw(ArgumentError("background diffusivities must be nonnegative"))
    relaxation_time_t > zero(T) || throw(ArgumentError("relaxation_time must be positive"))
    critical_ri_t > zero(T) || throw(ArgumentError("bulk_richardson_critical must be positive"))
    surface_mode in (:prescribed, :slab) || throw(ArgumentError(
        "surface_mode must be :prescribed or :slab",
    ))
    radiation_max_t >= zero(T) || throw(ArgumentError("radiation_max must be nonnegative"))
    longwave_loss_t >= zero(T) || throw(ArgumentError("longwave_loss must be nonnegative"))
    ground_conductance_t >= zero(T) || throw(ArgumentError("ground_conductance must be nonnegative"))
    return SCMConfig(dt_t, diff(collect(z)), z0_t, theta_flux_t, radiative_tendency_t,
        heat_capacity_t, km_t, kh_t, geostrophic_u_t, geostrophic_v_t, relaxation_time_t,
        critical_ri_t, surface_mode, radiation_max_t, longwave_loss_t, radiation_peak_time_t,
        soil_temperature_t, ground_conductance_t)
end

struct SCMDiagnostics{T<:AbstractFloat}
    t2m::T
    boundary_layer_height::T
    surface_sensible_heat_flux::T
    net_radiation::T
    ground_heat_flux::T
    gated_levels::Int
    lambda_f::Vector{T}
    gated_mask::BitVector
    ri_g::Vector{T}
    override_mask::BitVector
    gate_diagnostics::Vector{GateDiagnostics{T}}
end

"""
    net_radiation(config, elapsed_seconds)

Returns the periodic, clear-sky net radiation for the synthetic Cabauw surface.
The shortwave component peaks at `config.radiation_peak_time`; `longwave_loss`
provides the nocturnal radiative cooling term.
"""
function net_radiation(config::SCMConfig{T}, elapsed_seconds::T) where {T<:AbstractFloat}
    phase = T(2 * pi) * mod(elapsed_seconds - config.radiation_peak_time, T(86_400)) / T(86_400)
    return config.radiation_max * max(zero(T), cos(phase)) - config.longwave_loss
end

"""
    surface_sensible_heat_flux(state, K_h_surface)

Returns the upward sensible heat flux from the slab surface to the first model
level. Positive values warm the atmosphere and cool the surface.
"""
function surface_sensible_heat_flux(state::SCMState{T}, K_h_surface::T) where {T<:AbstractFloat}
    K_h_surface >= zero(T) || throw(ArgumentError("surface heat diffusivity must be nonnegative"))
    return T(RHO_CP_AIR) * K_h_surface *
        (state.skin_temperature - state.theta_v[1]) / state.z[1]
end

function initialize_cabauw_state(z::AbstractVector{T}) where {T<:AbstractFloat}
    z0 = T(0.15)
    u = T(4) .+ T(0.015) .* z
    v = T(0.5) .* exp.(-((z .- T(120)) ./ T(60)) .^ 2)
    theta_v = T(280) .+ T(0.015) .* z
    tke = T(0.08) .* exp.(-z ./ T(120))
    return SCMState(z, u, v, theta_v, tke, T(279))
end

function _vertical_gradient(values::Vector{T}, z::Vector{T}) where {T<:AbstractFloat}
    gradient = similar(values)
    gradient[1] = (values[2] - values[1]) / (z[2] - z[1])
    gradient[end] = (values[end] - values[end-1]) / (z[end] - z[end-1])
    for index in 2:(length(values)-1)
        gradient[index] = (values[index+1] - values[index-1]) / (z[index+1] - z[index-1])
    end
    return gradient
end

function solve_tridiagonal!(
    solution::Vector{T}, lower::Vector{T}, diagonal::Vector{T}, upper::Vector{T}, rhs::Vector{T},
) where {T<:AbstractFloat}
    n_levels = length(diagonal)
    length(solution) == n_levels && length(rhs) == n_levels || throw(ArgumentError(
        "solution, diagonal, and rhs must have matching lengths",
    ))
    length(lower) == n_levels - 1 && length(upper) == n_levels - 1 || throw(ArgumentError(
        "lower and upper diagonals must have length n - 1",
    ))
    for index in 2:n_levels
        pivot = diagonal[index-1]
        abs(pivot) > eps(T) || throw(ArgumentError("singular tridiagonal diffusion operator"))
        factor = lower[index-1] / pivot
        diagonal[index] -= factor * upper[index-1]
        rhs[index] -= factor * rhs[index-1]
    end
    abs(diagonal[end]) > eps(T) || throw(ArgumentError("singular tridiagonal diffusion operator"))
    solution[end] = rhs[end] / diagonal[end]
    for index in (n_levels-1):-1:1
        solution[index] = (rhs[index] - upper[index] * solution[index+1]) / diagonal[index]
    end
    return solution
end

function _implicit_diffuse!(
    field::Vector{T}, diffusivity::Vector{T}, z::Vector{T}, dt::T, lower_flux::T,
    lower_flux_gradient::T=zero(T),
) where {T<:AbstractFloat}
    n_levels = length(field)
    lower = zeros(T, n_levels - 1)
    diagonal = ones(T, n_levels)
    upper = zeros(T, n_levels - 1)
    rhs = copy(field)
    cell_widths = Vector{T}(undef, n_levels)
    cell_widths[1] = (z[2] - z[1]) / T(2)
    cell_widths[end] = (z[end] - z[end-1]) / T(2)
    for index in 2:(n_levels-1)
        cell_widths[index] = (z[index+1] - z[index-1]) / T(2)
    end
    for index in 1:(n_levels-1)
        spacing = z[index+1] - z[index]
        face_diffusivity = (diffusivity[index] + diffusivity[index+1]) / T(2)
        flux_coefficient = dt * face_diffusivity / spacing
        diagonal[index] += flux_coefficient / cell_widths[index]
        diagonal[index+1] += flux_coefficient / cell_widths[index+1]
        upper[index] -= flux_coefficient / cell_widths[index]
        lower[index] -= flux_coefficient / cell_widths[index+1]
    end
    rhs[1] += dt * lower_flux / cell_widths[1]
    diagonal[1] -= dt * lower_flux_gradient / cell_widths[1]
    return solve_tridiagonal!(field, lower, diagonal, upper, rhs)
end

function _boundary_layer_height(state::SCMState{T}, config::SCMConfig{T}) where {T<:AbstractFloat}
    theta_surface = state.theta_v[1]
    u_surface = state.u[1]
    v_surface = state.v[1]
    for index in 2:length(state.z)
        speed_difference_squared = (state.u[index] - u_surface)^2 + (state.v[index] - v_surface)^2 + eps(T)
        Ri_b = T(GRAVITY) / theta_surface * (state.theta_v[index] - theta_surface) *
               (state.z[index] - state.z[1]) / speed_difference_squared
        Ri_b >= config.bulk_richardson_critical && return state.z[index]
    end
    return state.z[end]
end

function _diagnose_t2m(state::SCMState{T}, config::SCMConfig{T}) where {T<:AbstractFloat}
    z1 = state.z[1]
    z1 > config.z0 || throw(ArgumentError("lowest model level must exceed z0"))
    interpolation_weight = log(T(2) / config.z0) / log(z1 / config.z0)
    return state.skin_temperature + interpolation_weight * (state.theta_v[1] - state.skin_temperature)
end

"""
    step_scm!(state, config, gating_params, gs_config;
              closure_mode=:full_gspt, enable_gating=true)

Advances the opt-in single-column benchmark one backward-Euler diffusion step.
`closure_mode=:control` applies classical Richardson-number quenching.
`:gate_only` preserves that closure but restores floors only for auditable gate
overrides. `:full_gspt` uses the regularized transport diagnostic and its
geometry-aware floor policy. All modes share the same implicit transport and
surface-energy solver.
"""
function step_scm!(
    state::SCMState{T},
    config::SCMConfig{T},
    gating_params::BifurcationGatingParams{T},
    gs_config::GSPTModelConfig{T};
    closure_mode::Symbol=:full_gspt,
    enable_gating::Bool=true,
) where {T<:AbstractFloat}
    closure_mode in (:control, :gate_only, :full_gspt) || throw(ArgumentError(
        "closure_mode must be :control, :gate_only, or :full_gspt",
    ))
    n_levels = length(state.z)
    length(config.dz) == n_levels - 1 || throw(ArgumentError("state and config grids differ"))
    shear = sqrt.(_vertical_gradient(state.u, state.z) .^ 2 + _vertical_gradient(state.v, state.z) .^ 2)
    N2 = T(GRAVITY) ./ state.theta_v .* _vertical_gradient(state.theta_v, state.z)
    zeta = N2 ./ (shear .^ 2 .+ sqrt(eps(T)))
    zeta_z = _vertical_gradient(zeta, state.z)
    base_K_m = gs_config.l0 .* sqrt.(max.(state.tke, zero(T)))
    base_K_h = base_K_m ./ T(1.5)
    K_m = similar(base_K_m)
    K_h = similar(base_K_h)
    lambda_f = Vector{T}(undef, n_levels)
    gated_mask = falses(n_levels)
    override_mask = falses(n_levels)
    gate_diagnostics = Vector{GateDiagnostics{T}}(undef, n_levels)

    for index in eachindex(state.z)
        jacobian = map_gabls3_to_jacobian(index, state.tke, shear, N2, gs_config)
        rig_shutdown = zeta[index] >= config.bulk_richardson_critical
        gate_diagnostics[index] = evaluate_gate_step!(
            state.gating_states[index], jacobian, zeta_z[index], rig_shutdown,
            gating_params; complex_policy=:real_part,
        )
        lambda_f[index] = gate_diagnostics[index].lambda_f
        override_mask[index] = gate_diagnostics[index].q_override
        stability_factor = rig_shutdown ? zero(T) :
            max(zero(T), one(T) - max(zeta[index], zero(T)) / config.bulk_richardson_critical)^2

        if closure_mode === :control || closure_mode === :gate_only
            K_m[index] = max(config.km_background, base_K_m[index] * stability_factor)
            K_h[index] = max(config.kh_background, base_K_h[index] * stability_factor)
            if closure_mode === :gate_only && enable_gating && override_mask[index]
                K_m[index] = max(K_m[index], gating_params.km_floor)
                K_h[index] = max(K_h[index], gating_params.kh_floor)
                gated_mask[index] = true
            end
        else
            K_m[index] = max(config.km_background, base_K_m[index])
            K_h[index] = max(config.kh_background, base_K_h[index])
            if enable_gating && gate_diagnostics[index].gate_active &&
               gate_diagnostics[index].is_coordinate_regular
                K_m[index] = max(K_m[index], gating_params.km_floor)
                K_h[index] = max(K_h[index], gating_params.kh_floor)
                gated_mask[index] = true
            end
        end
    end

    relaxation = config.dt / config.relaxation_time
    state.u .+= relaxation .* (config.geostrophic_u .- state.u)
    state.v .+= relaxation .* (config.geostrophic_v .- state.v)
    _implicit_diffuse!(state.u, K_m, state.z, config.dt, zero(T))
    _implicit_diffuse!(state.v, K_m, state.z, config.dt, zero(T))

    net_radiative_flux = zero(T)
    ground_heat_flux = zero(T)
    surface_heat_flux = if config.surface_mode === :slab
        net_radiative_flux = net_radiation(config, state.elapsed_seconds)
        heat_transfer = T(RHO_CP_AIR) * K_h[1] / state.z[1]
        surface_factor = config.dt / config.surface_heat_capacity
        surface_denominator = one(T) + surface_factor *
            (config.ground_conductance + heat_transfer)
        surface_constant = (state.skin_temperature + surface_factor *
            (net_radiative_flux + config.ground_conductance * config.soil_temperature)) /
            surface_denominator
        surface_theta_coefficient = surface_factor * heat_transfer / surface_denominator
        theta_flux_constant = K_h[1] * surface_constant / state.z[1]
        theta_flux_gradient = K_h[1] * (surface_theta_coefficient - one(T)) / state.z[1]
        _implicit_diffuse!(state.theta_v, K_h, state.z, config.dt, theta_flux_constant,
            theta_flux_gradient)
        state.skin_temperature = surface_constant + surface_theta_coefficient * state.theta_v[1]
        ground_heat_flux = config.ground_conductance *
            (state.skin_temperature - config.soil_temperature)
        surface_sensible_heat_flux(state, K_h[1])
    else
        T(RHO_CP_AIR) * config.theta_flux
    end
    if config.surface_mode === :prescribed
        theta_flux = surface_heat_flux / T(RHO_CP_AIR)
        _implicit_diffuse!(state.theta_v, K_h, state.z, config.dt, theta_flux)
    end

    if config.surface_mode === :slab
        nothing
    else
        state.skin_temperature += config.dt * (config.surface_radiative_tendency -
            surface_heat_flux / config.surface_heat_capacity)
    end
    state.elapsed_seconds += config.dt

    return SCMDiagnostics(
        _diagnose_t2m(state, config), _boundary_layer_height(state, config), surface_heat_flux,
        net_radiative_flux, ground_heat_flux, count(gated_mask), lambda_f, gated_mask,
        zeta, override_mask, gate_diagnostics,
    )
end

end # module GABLS3Stepper