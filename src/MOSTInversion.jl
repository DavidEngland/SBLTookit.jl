module MOSTInversion

export z0hr_split, compute_lagrange_coeffs, solve_zeta_z0hr

"""
    z0hr_split(Ri; ϵ=1e-6)

Zero-Offset Hyperbolic Regularization (Z0HR). Decomposes `Ri` into C^∞-differentiable
positive (stable) and negative (unstable) components without `if/else` branching.
"""
@inline function z0hr_split(Ri::T; ϵ::T=T(1e-6)) where {T<:AbstractFloat}
    # Hyperbolic regularized absolute value: |Ri|_ϵ = √(Ri² + ϵ²)
    abs_ri_ϵ = sqrt(Ri * Ri + ϵ * ϵ)

    # Smooth positive (stable) and negative (unstable) projections
    Ri_pos = T(0.5) * (Ri + abs_ri_ϵ - ϵ)
    Ri_neg = T(0.5) * (Ri - abs_ri_ϵ + ϵ)

    # Smooth C^∞ weight for regime blending: w ∈ (0, 1)
    w_stable = T(0.5) * (T(1.0) + Ri / abs_ri_ϵ)

    return Ri_pos, Ri_neg, w_stable
end

"""
    compute_lagrange_coeffs(Pr_t0, β_m, β_h)

Layer 1: Physics Engine. Derives forward Taylor coefficients c_k for Ri(ζ)
and computes universal Lagrange inversion coefficients a_k for ζ(Ri) near neutral (Ri ≈ 0).
"""
@inline function compute_lagrange_coeffs(Pr_t0::T, β_m::T, β_h::T) where {T<:AbstractFloat}
    # Forward profile expansion: Ri(ζ) = c₁ζ + c₂ζ² + c₃ζ³ + O(ζ⁴)
    c1 = Pr_t0
    c2 = β_h - T(2.0) * β_m * Pr_t0
    c3 = T(3.0) * β_m * β_m * Pr_t0 - T(2.0) * β_m * β_h

    # Layer 2: Universal Lagrange Inversion: ζ(Ri) = a₁Ri + a₂Ri² + a₃Ri³ + O(Ri⁴)
    c1_inv = T(1.0) / c1
    a1 = c1_inv
    a2 = -c2 * (c1_inv^3)
    a3 = (T(2.0) * c2 * c2 - c1 * c3) * (c1_inv^5)

    return (a1, a2, a3)
end

"""
    eval_lagrange_kernel(Ri, (a1, a2, a3))

Layer 2: Universal Lagrange Inversion Kernel using Horner's method for Ri ≈ 0.
"""
@inline function eval_lagrange_kernel(Ri::T, a::NTuple{3, T}) where {T<:AbstractFloat}
    a1, a2, a3 = a
    return Ri * (a1 + Ri * (a2 + Ri * a3))
end

"""
    solve_zeta_z0hr(Ri, a_coeffs; Ri_c=0.19, Pr_t0=1.0, ϵ=1e-6)

Unified GPU/SIMD branch-free solver evaluating ζ(Ri) across unstable, near-neutral,
and stable regimes via smooth Z0HR weighting.
"""
@inline function solve_zeta_z0hr(
    Ri::T,
    a_coeffs::NTuple{3, T};
    Ri_c::T=T(0.19),
    Pr_t0::T=T(1.0),
    ϵ::T=T(1e-6)
) where {T<:AbstractFloat}

    # 1. Branch-free Z0HR split
    Ri_pos, Ri_neg, w_stable = z0hr_split(Ri; ϵ=ϵ)

    # 2. Near-Neutral / Unstable Branch via Lagrange Series Reversion
    zeta_neutral = eval_lagrange_kernel(Ri, a_coeffs)

    # 3. Stable Branch via Quadratic Closure (S_m = (1 - Ri/Ri_c)²)
    # Regularize denominator to prevent division by zero near critical Ri_c
    denom_stable = max(T(1.0) - Ri_pos / Ri_c, ϵ)
    zeta_stable = (Ri_pos / Pr_t0) / denom_stable

    # 4. C^∞ smooth blend (No thread divergence on GPUs)
    zeta = (T(1.0) - w_stable) * zeta_neutral + w_stable * zeta_stable

    return zeta
end

end # module