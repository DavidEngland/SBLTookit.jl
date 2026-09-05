# =============================================================================
# SBLToolkit: src/sheba-bulk-fallback.jl
# Production-Grade Bulk Fallback Ingestion & Telemetry Processing Layer
# =============================================================================
# Tailored specifically for the SHEBA 2-level vertical tower footprint (2.5m, 10m).
# Automatically bypasses the differential geometric curvature operators (D2)
# to avoid stencil under-determination failure (StencilCollapseException).
#
# Enforces a mathematically rigorous Bulk Richardson (Ri_b) formulation with
# conditional regularized exchange floors, appending a STENCIL_DEGRADED telemetry
# flag to alert downstream analysis layers.
# =============================================================================

module SHEBABulkFallback

using LinearAlgebra
using Statistics
using Printf

export SHEBAParams, SHEBADIagnosticRow, evaluate_sheba_profile, run_sheba_climatology_batch

"""
    SHEBAParams{T<:AbstractFloat}

Configuration container for SHEBA bulk fallback diagnostics.
"""
struct SHEBAParams{T<:AbstractFloat}
    g::T            # Gravitational acceleration (m/s²)
    theta_ref::T    # Reference potential temperature (K)
    Ri_c::T         # Critical Richardson threshold for local collapse
    kh_floor::T     # Regularized exchange floor for heat (m²/s)
    km_floor::T     # Regularized exchange floor for momentum (m²/s)
end

function SHEBAParams(;
    g = 9.81,
    theta_ref = 255.0, # Typical cold polar SBL reference
    Ri_c = 0.25,
    kh_floor = 0.015,  # 0.015 m²/s target polar exchange floor
    km_floor = 0.10    # 0.10 m²/s target polar exchange floor
)
    return SHEBAParams{Float64}(g, theta_ref, Ri_c, kh_floor, km_floor)
end

"""
    SHEBADIagnosticRow{T<:AbstractFloat}

Represents the processed telemetry row of a SHEBA sounding.
"""
struct SHEBADIagnosticRow{T<:AbstractFloat}
    time_sec::T
    Ri_b::T
    K_h::T
    K_m::T
    STENCIL_DEGRADED::Bool
    gated_active::Bool
end

"""
    compute_bulk_richardson(z1, z2, u1, u2, v1, v2, theta1, theta2, params)

Computes the bulk Richardson number across two tower levels.
"""
function compute_bulk_richardson(
    z1::T, z2::T,
    u1::T, u2::T,
    v1::T, v2::T,
    theta1::T, theta2::T,
    params::SHEBAParams{T}
) where T<:AbstractFloat
    d_theta = theta2 - theta1
    du = u2 - u1
    dv = v2 - v1
    shear_sq = du^2 + dv^2
    
    # Regularization floor to prevent division-by-zero under absolute calm
    if shear_sq < 1e-6
        shear_sq = 1e-6
    end
    
    mean_theta = 0.5 * (theta1 + theta2)
    Ri_b = (params.g / mean_theta) * (z2 - z1) * d_theta / shear_sq
    return Ri_b
end

"""
    evaluate_sheba_profile(time, z, u, v, theta, params)

Processes a single 2-level SHEBA tower sounding. Resolves bulk exchange coefficients
and applies the regularized exchange floor conditionally under strong stability.
"""
function evaluate_sheba_profile(
    time::T,
    z::Vector{T},
    u::Vector{T},
    v::Vector{T},
    theta::Vector{T},
    params::SHEBAParams{T}
) where T<:AbstractFloat
    Nz = length(z)
    
    # 1. Enforce strict vertical grid size gate
    # If Nz < 3, differential geometric curvature cannot be resolved (STENCIL_DEGRADED = true)
    stencil_degraded = Nz < 3
    
    # 2. Extract bulk variables
    # For SHEBA, we map the two available levels: z1 = 2.5m, z2 = 10.0m
    z1, z2 = z[1], z[2]
    u1, u2 = u[1], u[2]
    v1, v2 = v[1], v[2]
    theta1, theta2 = theta[1], theta[2]
    
    # 3. Compute bulk Richardson stability
    Ri_b = compute_bulk_richardson(z1, z2, u1, u2, v1, v2, theta1, theta2, params)
    
    # 4. Standard unregularized local mixing calculation
    # Base mixing scales with wind shear and mixing length l0 = 0.40m
    l0 = 0.40
    S_bulk = sqrt((u2 - u1)^2 + (v2 - v1)^2) / (z2 - z1)
    K_m_raw = (l0^2) * S_bulk * max(0.0, 1.0 - Ri_b / params.Ri_c)^2
    K_h_raw = K_m_raw # Using Reynolds analogy for unregularized baseline
    
    # 5. Gating Decision: Apply polar exchange floors when Ri_b exceeds critical cutoff
    gated_active = Ri_b >= params.Ri_c
    
    K_m_eff = gated_active ? max(K_m_raw, params.km_floor) : K_m_raw
    K_h_eff = gated_active ? max(K_h_raw, params.kh_floor) : K_h_raw
    
    return SHEBADIagnosticRow{T}(time, Ri_b, K_h_eff, K_m_eff, stencil_degraded, gated_active)
end

"""
    run_sheba_climatology_batch(times, z_grid, u_data, v_data, theta_data, params)

Processes a multi-month SHEBA dataset, compiling bulk telemetry diagnostics.
"""
function run_sheba_climatology_batch(
    times::Vector{T},
    z_grid::Vector{T},
    u_data::Matrix{T},
    v_data::Matrix{T},
    theta_data::Matrix{T},
    params::SHEBAParams{T}
) where T<:AbstractFloat
    Nt = length(times)
    results = Vector{SHEBADIagnosticRow{T}}(undef, Nt)
    
    for t in 1:Nt
        results[t] = evaluate_sheba_profile(
            times[t], z_grid, u_data[:, t], v_data[:, t], theta_data[:, t], params
        )
    end
    
    return results
end

end # module SHEBABulkFallback
