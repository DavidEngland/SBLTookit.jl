module GABLS3Stepper

using LinearAlgebra
using ..SBLGating
using ..GABLS3Adapters

export SCMState, SCMConfig, SCMDiagnostics, solve_tridiagonal!, step_scm!,
    initialize_cabauw_state

const GRAVITY = 9.81

mutable struct SCMState{T<:AbstractFloat}
    z::Vector{T}
    u::Vector{T}
    v::Vector{T}
    theta_v::Vector{T}
    tke::Vector{T}
    skin_temperature::T
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
        skin_temperature, [GatingState(T) for _ in 1:n_levels],
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
) where {T<:AbstractFloat}
    length(z) >= 2 && issorted(z) && all(diff(z) .> zero(T)) || throw(ArgumentError(
        "SCM configuration requires strictly increasing levels",
    ))
    values = T.((dt, z0, theta_flux, surface_radiative_tendency, surface_heat_capacity,
        km_background, kh_background, geostrophic_u, geostrophic_v, relaxation_time,
        bulk_richardson_critical))
    all(isfinite, values) || throw(ArgumentError("SCM configuration values must be finite"))
    dt_t, z0_t, theta_flux_t, radiative_tendency_t, heat_capacity_t, km_t, kh_t,
        geostrophic_u_t, geostrophic_v_t, relaxation_time_t, critical_ri_t = values
    dt_t > zero(T) || throw(ArgumentError("dt must be positive"))
    z0_t > zero(T) || throw(ArgumentError("z0 must be positive"))
    heat_capacity_t > zero(T) || throw(ArgumentError("surface_heat_capacity must be positive"))
    km_t >= zero(T) && kh_t >= zero(T) || throw(ArgumentError("background diffusivities must be nonnegative"))
    relaxation_time_t > zero(T) || throw(ArgumentError("relaxation_time must be positive"))
    critical_ri_t > zero(T) || throw(ArgumentError("bulk_richardson_critical must be positive"))
    return SCMConfig(dt_t, diff(collect(z)), z0_t, theta_flux_t, radiative_tendency_t,
        heat_capacity_t, km_t, kh_t, geostrophic_u_t, geostrophic_v_t, relaxation_time_t,
        critical_ri_t)
end

struct SCMDiagnostics{T<:AbstractFloat}
    t2m::T
    boundary_layer_height::T
    surface_sensible_heat_flux::T
    gated_levels::Int
    lambda_f::Vector{T}
    gated_mask::BitVector
end

function initialize_cabauw_state(z::AbstractVector{T}) where {T<:AbstractFloat}
    z0 = T(0.15)
    u = T(4) .+ T(0.015) .* z
    v = T(0.5) .* exp.(-((z .- T(120)) ./ T(60)).^2)
    theta_v = T(280) .+ T(0.015) .* z
    tke = T(0.08) .* exp.(-z ./ T(120))
    return SCMState(z, u, v, theta_v, tke, T(279))
end

function _vertical_gradient(values::Vector{T}, z::Vector{T}) where {T<:AbstractFloat}
    gradient = similar(values)
    gradient[1] = (values[2] - values[1]) / (z[2] - z[1])
    gradient[end] = (values[end] - values[end - 1]) / (z[end] - z[end - 1])
    for index in 2:(length(values) - 1)
        gradient[index] = (values[index + 1] - values[index - 1]) / (z[index + 1] - z[index - 1])
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
        pivot = diagonal[index - 1]
        abs(pivot) > eps(T) || throw(ArgumentError("singular tridiagonal diffusion operator"))
        factor = lower[index - 1] / pivot
        diagonal[index] -= factor * upper[index - 1]
        rhs[index] -= factor * rhs[index - 1]
    end
    abs(diagonal[end]) > eps(T) || throw(ArgumentError("singular tridiagonal diffusion operator"))
    solution[end] = rhs[end] / diagonal[end]
    for index in (n_levels - 1):-1:1
        solution[index] = (rhs[index] - upper[index] * solution[index + 1]) / diagonal[index]
    end
    return solution
end

function _implicit_diffuse!(
    field::Vector{T}, diffusivity::Vector{T}, z::Vector{T}, dt::T, lower_flux::T,
) where {T<:AbstractFloat}
    n_levels = length(field)
    lower = zeros(T, n_levels - 1)
    diagonal = ones(T, n_levels)
    upper = zeros(T, n_levels - 1)
    rhs = copy(field)
    for index in 1:(n_levels - 1)
        spacing = z[index + 1] - z[index]
        face_diffusivity = (diffusivity[index] + diffusivity[index + 1]) / T(2)
        coefficient = dt * face_diffusivity / spacing
        if index == 1
            upper[index] -= coefficient / spacing
            diagonal[index] += coefficient / spacing
            rhs[index] += dt * lower_flux / spacing
        else
            lower[index - 1] -= coefficient / spacing
            diagonal[index] += coefficient / spacing
            upper[index] -= coefficient / spacing
            diagonal[index + 1] += coefficient / spacing
            lower[index] -= coefficient / spacing
        end
    end
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
    step_scm!(state, config, gating_params, gs_config; enable_gating=true)

Advances the opt-in single-column benchmark one backward-Euler diffusion step.
The gate only replaces diagnosed diffusivities after they are computed; it does
not change the underlying stability or turbulence closure.
"""
function step_scm!(
    state::SCMState{T},
    config::SCMConfig{T},
    gating_params::BifurcationGatingParams{T},
    gs_config::GSPTModelConfig{T};
    enable_gating::Bool=true,
) where {T<:AbstractFloat}
    n_levels = length(state.z)
    length(config.dz) == n_levels - 1 || throw(ArgumentError("state and config grids differ"))
    shear = sqrt.(_vertical_gradient(state.u, state.z).^2 + _vertical_gradient(state.v, state.z).^2)
    N2 = T(GRAVITY) ./ state.theta_v .* _vertical_gradient(state.theta_v, state.z)
    zeta = N2 ./ (shear.^2 .+ sqrt(eps(T)))
    zeta_z = _vertical_gradient(zeta, state.z)
    K_m = max.(config.km_background, gs_config.l0 .* sqrt.(max.(state.tke, zero(T))))
    K_h = max.(config.kh_background, K_m ./ T(1.5))
    lambda_f = Vector{T}(undef, n_levels)
    gated_mask = falses(n_levels)

    for index in eachindex(state.z)
        jacobian = map_gabls3_to_jacobian(index, state.tke, shear, N2, gs_config)
        lambda_f[index] = extract_fast_eigenvalue(jacobian, state.gating_states[index])
        if enable_gating
            K_m[index], K_h[index], gated_mask[index] = update_gating_state!(
                state.gating_states[index], lambda_f[index], zeta_z[index], K_m[index], K_h[index], gating_params,
            )
        end
    end

    relaxation = config.dt / config.relaxation_time
    state.u .+= relaxation .* (config.geostrophic_u .- state.u)
    state.v .+= relaxation .* (config.geostrophic_v .- state.v)
    _implicit_diffuse!(state.u, K_m, state.z, config.dt, zero(T))
    _implicit_diffuse!(state.v, K_m, state.z, config.dt, zero(T))
    _implicit_diffuse!(state.theta_v, K_h, state.z, config.dt, config.theta_flux)
    surface_heat_flux = T(1.225 * 1004.67) * config.theta_flux
    state.skin_temperature += config.dt * (config.surface_radiative_tendency -
        surface_heat_flux / config.surface_heat_capacity)

    return SCMDiagnostics(
        _diagnose_t2m(state, config), _boundary_layer_height(state, config), surface_heat_flux,
        count(gated_mask), lambda_f, gated_mask,
    )
end

end # module GABLS3Stepper