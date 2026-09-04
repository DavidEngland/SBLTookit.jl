module SBLGating

using LinearAlgebra

export BifurcationGatingParams, GatingState, extract_fast_eigenvalue,
    update_gating_state!

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
end

function BifurcationGatingParams(
    ; epsilon_on::Real=-0.10,
    epsilon_off::Real=-0.05,
    zeta_z_tol::Real=1e-2,
    km_floor::Real=0.1,
    kh_floor::Real=0.01,
)
    T = promote_type(typeof(float(epsilon_on)), typeof(float(epsilon_off)),
        typeof(float(zeta_z_tol)), typeof(float(km_floor)), typeof(float(kh_floor)))
    values = T.((epsilon_on, epsilon_off, zeta_z_tol, km_floor, kh_floor))
    all(isfinite, values) || throw(ArgumentError("gating parameters must be finite"))
    epsilon_on_t, epsilon_off_t, zeta_z_tol_t, km_floor_t, kh_floor_t = values
    epsilon_on_t < epsilon_off_t < zero(T) || throw(ArgumentError(
        "require epsilon_on < epsilon_off < 0 for attracting-mode hysteresis",
    ))
    zeta_z_tol_t >= zero(T) || throw(ArgumentError("zeta_z_tol must be nonnegative"))
    km_floor_t >= zero(T) || throw(ArgumentError("km_floor must be nonnegative"))
    kh_floor_t >= zero(T) || throw(ArgumentError("kh_floor must be nonnegative"))
    return BifurcationGatingParams(
        epsilon_on_t,
        epsilon_off_t,
        zeta_z_tol_t,
        km_floor_t,
        kh_floor_t,
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