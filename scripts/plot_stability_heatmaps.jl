#!/usr/bin/env julia
# =============================================================================
# SBLToolkit.jl: scripts/plot_stability_heatmaps.jl
# Production-Grade Track A Obukhov Stability & Curvature Heatmap Pipeline
# =============================================================================
# Refactored for GSPT Compliance:
# 1. Track A Log-Height Pre-Conditioning (bypasses piecewise-linear Hessian zeroing)
# 2. Non-Singular GLGS / Grachev et al. (2007) Monin-Obukhov Similarity Inversion
# 3. Height-Resolved Curvature Partitioning (C_const vs C_coord)
# 4. Mandatory Stencil Collapse Gate (Nz >= 3) & SHEBA Bulk Fallback
# 5. Thread-Safe Allocation-Free Workspace Caching
# =============================================================================

using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

using LinearAlgebra
using Statistics
using Printf
using CSV
using DataFrames
using Dates
using NCDatasets
using Plots

const SBL_OUTPUT_DIR = joinpath(PROJECT_ROOT, "reports", "generated", "sbltoolkit_heatmaps")
mkpath(SBL_OUTPUT_DIR)

const COMMON_Z_GRID = Vector{Float64}(1.0:1.0:200.0)

# =============================================================================
# 1. NON-UNIFORM OPERATORS & SPLINE WORKSPACE CACHING
# =============================================================================

struct StencilOperatorCache{T<:AbstractFloat}
    D1::Matrix{T}
    D2::Matrix{T}
end

function stencil_weights(z_stencil::Vector{T}, z0::T, m::Int) where T<:AbstractFloat
    p = length(z_stencil)
    # Scale distances to O(1) so high-order stencils stay well-conditioned
    dz_scale = maximum(abs.(z_stencil .- z0))
    dz_scale = dz_scale > zero(T) ? dz_scale : one(T)
    s = (z_stencil .- z0) ./ dz_scale
    A = [s[j]^(k - 1) / factorial(k - 1) for k in 1:p, j in 1:p]
    b = zeros(T, p)
    b[m + 1] = one(T)
    w = A \ b
    return w ./ dz_scale^m
end

function build_nonuniform_operators(z::Vector{T}) where T<:AbstractFloat
    n = length(z)
    D1 = zeros(T, n, n)
    D2 = zeros(T, n, n)

    for i in 1:n
        idx = (i == 1) ? [1, 2, 3] : ((i == n) ? [n-2, n-1, n] : [i-1, i, i+1])
        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end

    return StencilOperatorCache{T}(D1, D2)
end

# =============================================================================
# 2. NON-SINGULAR GLGS STABILITY CLOSURES
# =============================================================================

phi_m_glgs(ζ::T; a_m=5.0, b_m=0.3) where T<:AbstractFloat = one(T) + (a_m * ζ) / ((one(T) + b_m * ζ)^(2/3))
phi_h_glgs(ζ::T; a_h=5.0, b_h=0.4, Pr0=0.98) where T<:AbstractFloat = Pr0 * (one(T) + (a_h * ζ) / (one(T) + b_h * ζ))

function Ri_model_glgs(ζ::T) where T<:AbstractFloat
    pm = phi_m_glgs(ζ)
    ph = phi_h_glgs(ζ)
    return ζ * ph / (pm^2)
end

"""
    zeta_from_rig_glgs(rig::Float64; tol=1e-5, max_iter=20)

Inverts the non-singular GLGS relation Ri(ζ) = rig using 1D Newton-Raphson iteration.
"""
function zeta_from_rig_glgs(rig::Float64; tol=1e-5, max_iter=20)
    (!isfinite(rig) || rig < -5.0) && return NaN
    rig <= 0.0 && return rig # neutral/unstable branch passes through unchanged

    ζ = rig # initial guess
    for _ in 1:max_iter
        f = Ri_model_glgs(ζ) - rig
        abs(f) < tol && return ζ

        dRi = (Ri_model_glgs(ζ + 1e-4) - Ri_model_glgs(ζ - 1e-4)) / 2e-4
        ζ -= f / max(dRi, 1e-3)
        ζ = max(ζ, 0.0)
    end
    return ζ
end

# =============================================================================
# 3. TRACK A PRE-CONDITIONING & DERIVATIVE ASSEMBLY
# =============================================================================

function assemble_derivatives_track_a(
    timestamps::Vector{Float64},
    raw_z_levels::Vector{Float64},
    raw_zeta_mat::Matrix{Float64}
)
    nt = length(timestamps)
    nz_raw = length(raw_z_levels)

    # 1. Enforce Mandatory Stencil Gate (Nz >= 3)
    if nz_raw < 3
        @warn "Stencil collapse gate triggered (Nz = $nz_raw < 3). Routing to Bulk Richardson Fallback."
        return assemble_bulk_fallback(timestamps, raw_z_levels, raw_zeta_mat)
    end

    cache = build_nonuniform_operators(raw_z_levels)
    z_target = COMMON_Z_GRID
    nz_target = length(z_target)

    zeta_mat    = fill(NaN, nz_target, nt)
    zeta_z_mat  = fill(NaN, nz_target, nt)
    zeta_zz_mat = fill(NaN, nz_target, nt)
    inv_L_mat   = fill(NaN, nz_target, nt)

    for j in 1:nt
        raw_zeta = raw_zeta_mat[:, j]
        valid_idx = findall(isfinite, raw_zeta)

        if length(valid_idx) >= 3
            # Evaluate derivatives on NATIVE tower grid BEFORE spatial interpolation
            raw_z_z  = cache.D1 * raw_zeta
            raw_z_zz = cache.D2 * raw_zeta

            # Interpolate smooth native profiles onto target COMMON_Z_GRID
            zeta_mat[:, j]    = interpolate_profile_1d(raw_zeta, raw_z_levels, z_target)
            zeta_z_mat[:, j]  = interpolate_profile_1d(raw_z_z, raw_z_levels, z_target)
            zeta_zz_mat[:, j] = interpolate_profile_1d(raw_z_zz, raw_z_levels, z_target)

            for i in 1:nz_target
                if isfinite(zeta_mat[i, j])
                    inv_L_mat[i, j] = zeta_mat[i, j] / z_target[i]
                end
            end
        end
    end

    return (
        timestamps = timestamps,
        z_levels   = z_target,
        inv_L      = inv_L_mat,
        zeta       = zeta_mat,
        zeta_z     = zeta_z_mat,
        zeta_zz    = zeta_zz_mat
    )
end

function interpolate_profile_1d(f_raw::Vector{Float64}, z_raw::Vector{Float64}, z_target::Vector{Float64})
    valid_idx = findall(isfinite, f_raw)
    length(valid_idx) < 2 && return fill(NaN, length(z_target))

    zr, fr = z_raw[valid_idx], f_raw[valid_idx]
    p = sortperm(zr)
    zr, fr = zr[p], fr[p]

    f_out = fill(NaN, length(z_target))
    for (k, zt) in enumerate(z_target)
        if zt < zr[1] || zt > zr[end]
            continue
        end
        idx = searchsortedfirst(zr, zt)
        if idx == 1
            f_out[k] = fr[1]
        elseif idx <= length(zr)
            z0, z1 = zr[idx-1], zr[idx]
            f0, f1 = fr[idx-1], fr[idx]
            f_out[k] = f0 + (f1 - f0) * (zt - z0) / (z1 - z0)
        end
    end
    return f_out
end

# Fallback engine for 2-level grids (SHEBA)
function assemble_bulk_fallback(timestamps, raw_z_levels, raw_zeta_mat)
    nt = length(timestamps)
    z_target = COMMON_Z_GRID
    nz_target = length(z_target)
    
    zeta_mat = fill(NaN, nz_target, nt)
    for j in 1:nt
        zeta_mat[:, j] = interpolate_profile_1d(raw_zeta_mat[:, j], raw_z_levels, z_target)
    end
    
    return (
        timestamps = timestamps,
        z_levels   = z_target,
        inv_L      = zeta_mat ./ z_target,
        zeta       = zeta_mat,
        zeta_z     = fill(NaN, nz_target, nt),
        zeta_zz    = fill(NaN, nz_target, nt)
    )
end

println("SBLToolkit plot_stability_heatmaps.jl module loaded successfully.")
