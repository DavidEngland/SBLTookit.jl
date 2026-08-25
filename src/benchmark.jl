#!/usr/bin/env julia

# 1. Define Module First
module MOSTInversion

export z0hr_split, compute_lagrange_coeffs, solve_zeta_z0hr

"""
    z0hr_split(Ri; ϵ=1e-6)

Zero-Offset Hyperbolic Regularization (Z0HR). Decomposes `Ri` into C^∞-differentiable
positive (stable) and negative (unstable) components without `if/else` branching.
"""
@inline function z0hr_split(Ri::T; ϵ::T=T(1e-6)) where {T<:AbstractFloat}
    abs_ri_ϵ = sqrt(Ri * Ri + ϵ * ϵ)
    Ri_pos = T(0.5) * (Ri + abs_ri_ϵ - ϵ)
    Ri_neg = T(0.5) * (Ri - abs_ri_ϵ + ϵ)
    w_stable = T(0.5) * (T(1.0) + Ri / abs_ri_ϵ)
    return Ri_pos, Ri_neg, w_stable
end

"""
    compute_lagrange_coeffs(Pr_t0, β_m, β_h)

Layer 1: Physics Engine. Derives forward Taylor coefficients c_k for Ri(ζ)
and computes universal Lagrange inversion coefficients a_k for ζ(Ri) near neutral (Ri ≈ 0).
"""
@inline function compute_lagrange_coeffs(Pr_t0::T, β_m::T, β_h::T) where {T<:AbstractFloat}
    c1 = Pr_t0
    c2 = β_h - T(2.0) * β_m * Pr_t0
    c3 = T(3.0) * β_m * β_m * Pr_t0 - T(2.0) * β_m * β_h

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
@inline function eval_lagrange_kernel(Ri::T, a::NTuple{3,T}) where {T<:AbstractFloat}
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
    a_coeffs::NTuple{3,T};
    Ri_c::T=T(0.19),
    Pr_t0::T=T(1.0),
    ϵ::T=T(1e-6)
) where {T<:AbstractFloat}

    Ri_pos, Ri_neg, w_stable = z0hr_split(Ri; ϵ=ϵ)
    zeta_neutral = eval_lagrange_kernel(Ri, a_coeffs)

    denom_stable = max(T(1.0) - Ri_pos / Ri_c, ϵ)
    zeta_stable = (Ri_pos / Pr_t0) / denom_stable

    zeta = (T(1.0) - w_stable) * zeta_neutral + w_stable * zeta_stable
    return zeta
end

end # module MOSTInversion

# 2. Import Module & Run Driver
using .MOSTInversion
using BenchmarkTools

Pr_t0 = 1.0
β_m = 4.7
β_h = 4.7

a_coeffs = compute_lagrange_coeffs(Pr_t0, β_m, β_h)

Ri_profile = Float64[-0.5, -0.05, 0.0, 0.02, 0.1, 0.18]
zeta_profile = zeros(Float64, length(Ri_profile))

# Explicit scalar reference broadcast prevents `@.` macro expansion bugs
zeta_profile .= solve_zeta_z0hr.(Ri_profile, Ref(a_coeffs))

println("Gradient Richardson Numbers (Ri):")
println(Ri_profile)
println("\nDerived Stability Parameters (ζ = z/L):")
println(zeta_profile)

println("\nPerformance Benchmark:")
@btime $zeta_profile .= solve_zeta_z0hr.($Ri_profile, Ref($a_coeffs))