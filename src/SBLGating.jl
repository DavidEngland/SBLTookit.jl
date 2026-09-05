module SBLGating

using LinearAlgebra

export BifurcationGatingParams, GatingState, GateDiagnostics,
    extract_fast_eigenvalue, evaluate_gate_step!, update_gating_state!

"""
    BifurcationGatingParams{T}

Parameters for a geometry-aware mixing-floor experiment. `epsilon_on` and
`epsilon_off` are negative fast-eigenvalue thresholds with
`epsilon_on < epsilon_off < 0`: the dynamical gate activates at or below the
former and releases at or above the latter.
"""
struct BifurcationGatingParams{T<:AbstractFloat}
    epsilon_on::T
    epsilon_off::T
    zeta_z_tol::T
    km_floor::T
    kh_floor::T
    min_mode_overlap::T
end

function BifurcationGatingParams(
    ; epsilon_on::Real=-0.10,
    epsilon_off::Real=-0.05,
    zeta_z_tol::Real=1e-2,
    km_floor::Real=0.1,
    kh_floor::Real=0.01,
    min_mode_overlap::Real=0.70,
)
    T = promote_type(typeof(float(epsilon_on)), typeof(float(epsilon_off)),
        typeof(float(zeta_z_tol)), typeof(float(km_floor)), typeof(float(kh_floor)),
        typeof(float(min_mode_overlap)))
    values = T.((epsilon_on, epsilon_off, zeta_z_tol, km_floor, kh_floor, min_mode_overlap))
    all(isfinite, values) || throw(ArgumentError("gating parameters must be finite"))
    epsilon_on_t, epsilon_off_t, zeta_z_tol_t, km_floor_t, kh_floor_t, min_mode_overlap_t = values
    epsilon_on_t < epsilon_off_t < zero(T) || throw(ArgumentError(
        "require epsilon_on < epsilon_off < 0 for attracting-mode hysteresis",
    ))
    zeta_z_tol_t >= zero(T) || throw(ArgumentError("zeta_z_tol must be nonnegative"))
    km_floor_t >= zero(T) || throw(ArgumentError("km_floor must be nonnegative"))
    kh_floor_t >= zero(T) || throw(ArgumentError("kh_floor must be nonnegative"))
    zero(T) <= min_mode_overlap_t <= one(T) || throw(ArgumentError(
        "min_mode_overlap must lie in [0, 1]",
    ))
    return BifurcationGatingParams(
        epsilon_on_t,
        epsilon_off_t,
        zeta_z_tol_t,
        km_floor_t,
        kh_floor_t,
        min_mode_overlap_t,
    )
end

"""
    GatingState(T)

Persistent state for hysteresis and fast-mode continuation at one vertical level.
"""
mutable struct GatingState{T<:AbstractFloat}
    is_active::Bool
    prev_v_f::Union{Nothing, Vector{T}}
end

GatingState(::Type{T}) where {T<:AbstractFloat} = GatingState{T}(false, nothing)

"""
    GateDiagnostics{T}

Immutable audit record for one gate evaluation. `q_override` is true only when
a Richardson-number shutdown request coincides with a continuous, active fast
mode at a locally regular coordinate mapping.
"""
struct GateDiagnostics{T<:AbstractFloat}
    lambda_f::T
    mode_overlap::T
    gate_active::Bool
    is_coordinate_regular::Bool
    q_override::Bool
end

function _normalized_real_eigensystem(J::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(J) == (2, 2) || throw(ArgumentError("fast-mode extraction requires a 2x2 Jacobian"))
    all(isfinite, J) || throw(ArgumentError("Jacobian entries must be finite"))
    decomposition = eigen(Matrix(J))
    all(isreal, decomposition.values) || throw(ArgumentError(
        "fast-mode extraction requires real eigenvalues; complex spectra need a declared policy",
    ))
    values = T.(real.(decomposition.values))
    vectors = T.(real.(decomposition.vectors))
    all(isfinite, vectors) || throw(ArgumentError("Jacobian eigenvectors must be finite"))
    for column in axes(vectors, 2)
        magnitude = norm(view(vectors, :, column))
        magnitude > zero(T) || throw(ArgumentError("Jacobian has a zero eigenvector"))
        vectors[:, column] ./= magnitude
    end
    return values, vectors
end

"""
    extract_fast_eigenvalue(J, state; complex_policy=:reject)

Returns the real eigenvalue of the tracked fast branch. On the first call, the
branch with the shortest timescale (largest absolute real eigenvalue) initializes
the track. Later calls select the eigenvector with the greatest absolute overlap
with the prior branch, avoiding ordering changes at crossings. Complex spectra
have no unique real fast eigenvector: `:reject` throws, while `:real_part`
returns the complex branch decay rate and clears vector continuation.
"""
function extract_fast_eigenvalue(
    J::AbstractMatrix{T},
    state::GatingState{T},
    ; complex_policy::Symbol=:reject,
) where {T<:AbstractFloat}
    complex_policy in (:reject, :real_part) || throw(ArgumentError(
        "complex_policy must be :reject or :real_part",
    ))
    size(J) == (2, 2) || throw(ArgumentError("fast-mode extraction requires a 2x2 Jacobian"))
    all(isfinite, J) || throw(ArgumentError("Jacobian entries must be finite"))
    decomposition = eigen(Matrix(J))
    if !all(isreal, decomposition.values)
        complex_policy === :reject && throw(ArgumentError(
            "fast-mode extraction requires real eigenvalues; complex spectra need a declared policy",
        ))
        real_values = T.(real.(decomposition.values))
        state.prev_v_f = nothing
        return real_values[argmax(abs.(real_values))]
    end
    values, vectors = _normalized_real_eigensystem(J)
    index = if isnothing(state.prev_v_f)
        argmax(abs.(values))
    else
        overlaps = abs.(transpose(state.prev_v_f) * vectors)
        argmax(vec(overlaps))
    end
    state.prev_v_f = copy(vectors[:, index])
    return values[index]
end

"""
    evaluate_gate_step!(state, J, zeta_z, rig_shutdown, params; complex_policy=:reject)

Evaluates the tracked fast mode, its continuity with the preceding eigenvector,
the attracting-mode hysteresis state, coordinate regularity, and the auditable
Richardson-shutdown override condition.
"""
function evaluate_gate_step!(
    state::GatingState{T},
    J::AbstractMatrix{T},
    zeta_z::T,
    rig_shutdown::Bool,
    params::BifurcationGatingParams{T};
    complex_policy::Symbol=:reject,
) where {T<:AbstractFloat}
    isfinite(zeta_z) || throw(ArgumentError("zeta_z must be finite"))
    previous_vector = state.prev_v_f
    lambda_f = extract_fast_eigenvalue(J, state; complex_policy)
    mode_overlap = isnothing(previous_vector) || isnothing(state.prev_v_f) ? one(T) :
        abs(dot(previous_vector, state.prev_v_f))
    if mode_overlap < params.min_mode_overlap
        state.is_active = false
    elseif lambda_f <= params.epsilon_on
        state.is_active = true
    elseif lambda_f >= params.epsilon_off
        state.is_active = false
    end
    is_coordinate_regular = abs(zeta_z) <= params.zeta_z_tol
    q_override = rig_shutdown && state.is_active && is_coordinate_regular
    return GateDiagnostics(lambda_f, mode_overlap, state.is_active, is_coordinate_regular, q_override)
end

"""
    update_gating_state!(state, lambda_f, zeta_z, K_m, K_h, params)

Updates the attracting-mode hysteresis state and applies mixing floors only if
the state is active and the coordinate mapping is locally regular. This is a
numerical test of the geometry-dynamics decoupling hypothesis.
"""
function update_gating_state!(
    state::GatingState{T},
    lambda_f::T,
    zeta_z::T,
    K_m::T,
    K_h::T,
    params::BifurcationGatingParams{T},
) where {T<:AbstractFloat}
    all(isfinite, (lambda_f, zeta_z, K_m, K_h)) || throw(ArgumentError(
        "gating inputs must be finite",
    ))
    if lambda_f <= params.epsilon_on
        state.is_active = true
    elseif lambda_f >= params.epsilon_off
        state.is_active = false
    end

    is_gated = state.is_active && abs(zeta_z) <= params.zeta_z_tol
    if is_gated
        return max(K_m, params.km_floor), max(K_h, params.kh_floor), true
    end
    return K_m, K_h, false
end

end # module SBLGating