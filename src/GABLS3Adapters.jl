module GABLS3Adapters

export GSPTModelConfig, map_gabls3_to_jacobian, compute_fold_ratio

"""
    GSPTModelConfig{T}

Configuration for the reduced GABLS3 fast-slow diagnostic Jacobian.
"""
struct GSPTModelConfig{T<:AbstractFloat}
    epsilon::T
    epsilon_c::T
    l0::T
end

function GSPTModelConfig(
    ; epsilon::Real=1e-3,
    epsilon_c::Real=1e-3,
    l0::Real=15.0,
)
    T = promote_type(typeof(float(epsilon)), typeof(float(epsilon_c)), typeof(float(l0)))
    epsilon_t, epsilon_c_t, l0_t = T.((epsilon, epsilon_c, l0))
    all(isfinite, (epsilon_t, epsilon_c_t, l0_t)) || throw(ArgumentError(
        "GSPT model parameters must be finite",
    ))
    epsilon_t > zero(T) || throw(ArgumentError("epsilon must be positive"))
    epsilon_c_t > zero(T) || throw(ArgumentError("epsilon_c must be positive"))
    l0_t > zero(T) || throw(ArgumentError("l0 must be positive"))
    return GSPTModelConfig(epsilon_t, epsilon_c_t, l0_t)
end

"""
    map_gabls3_to_jacobian(z_idx, E, S, N2, config)

Constructs the scaled two-variable Jacobian
`[(dF/de)/epsilon (dF/dS)/epsilon; dg/de dg/dS]` at `z_idx`. The buoyancy
term uses the C1 Hill-type regularization proposed for the laminar limit.
"""
function map_gabls3_to_jacobian(
    z_idx::Integer,
    E::AbstractVector{T},
    S::AbstractVector{T},
    N2::AbstractVector{T},
    config::GSPTModelConfig{T},
) where {T<:AbstractFloat}
    n_levels = length(E)
    length(S) == n_levels && length(N2) == n_levels || throw(ArgumentError(
        "E, S, and N2 must have the same length",
    ))
    checkbounds(E, z_idx)
    e_raw, s_local, n2_local = E[z_idx], S[z_idx], N2[z_idx]
    all(isfinite, (e_raw, s_local, n2_local)) || throw(ArgumentError(
        "profile values must be finite",
    ))

    e_local = max(e_raw, sqrt(eps(T)))
    delta_reg = config.epsilon_c
    buoyancy_cap = config.l0 * n2_local
    d_buoyancy_de = buoyancy_cap * (2 * e_local * delta_reg^2) /
        (e_local^2 + delta_reg^2)^2
    d_dissipation_de = (T(1.5) * sqrt(e_local)) / config.l0

    F_e = config.l0 * s_local^2 - d_buoyancy_de - d_dissipation_de
    F_s = 2 * config.l0 * e_local * s_local
    gamma_s = T(0.1)
    r_s = T(0.01)
    g_e = -gamma_s * s_local
    g_s = -gamma_s * e_local - r_s
    inverse_epsilon = one(T) / config.epsilon

    return T[F_e * inverse_epsilon F_s * inverse_epsilon; g_e g_s]
end

"""
    compute_fold_ratio(zeta_z, zeta_zz, C_constitutive)

Returns a bounded mapping-curvature fraction in `[0, 1]`; values near one
indicate mapping-dominated curvature and are retained as an audit diagnostic.
"""
function compute_fold_ratio(
    zeta_z::T,
    zeta_zz::T,
    C_constitutive::T,
) where {T<:AbstractFloat}
    all(isfinite, (zeta_z, zeta_zz, C_constitutive)) || throw(ArgumentError(
        "fold-ratio inputs must be finite",
    ))
    mapping_curvature = abs(zeta_zz) / (abs(zeta_z) + sqrt(eps(T)))
    return mapping_curvature / (mapping_curvature + abs(C_constitutive) + eps(T))
end

end # module GABLS3Adapters