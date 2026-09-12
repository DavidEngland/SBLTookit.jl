#!/usr/bin/env julia
# =============================================================================
# C-Z0HR BENCHMARK v3: IN-PLACE ZERO-ALLOCATION MUTATING OPERATOR & C^∞ CLOSURE
# =============================================================================
# Demonstrates in-place, zero-allocation C-Z0HR mapping operator (compute_cz0hr_mapping!)
# paired with end-to-end C^∞ hyperbolic soft-plus closure in a 1D Low-Level Jet (LLJ).
# =============================================================================

using Printf
using LinearAlgebra

# Try to use Plots.jl for visualization if available
try
    using Plots
catch e
    @warn "Plots.jl package not found in current environment. Script will execute numerical benchmark only."
end

"""
    CZ0HRConfig{T<:AbstractFloat}

Immutable configuration container for C-Z0HR regularization parameters.
Type-parameterized on T to guarantee zero runtime boxing or dynamic dispatch.
"""
struct CZ0HRConfig{T<:AbstractFloat}
    Ri_c::T       # Critical Richardson number (default: 0.20)
    alpha::T      # Damping control parameter (default: 2.0)
    epsilon_c::T  # Softening term cap factor (default: 0.05)
    eps_hyper::T  # Hyperbolic C^∞ smoothing parameter (default: 1e-3)
end

function CZ0HRConfig(;
    Ri_c::T = 0.20,
    alpha::T = 2.0,
    epsilon_c::T = 0.05,
    eps_hyper::T = 1e-3
) where {T<:AbstractFloat}
    return CZ0HRConfig{T}(Ri_c, alpha, epsilon_c, eps_hyper)
end

"""
    Phi_eps(x::T, eps_hyper::T)

Smooth C^∞ hyperbolic distance metric replacing non-differentiable |x| kinks.
Annotated with @inline for LLVM vector register optimization.
"""
@inline Phi_eps(x::T, eps_hyper::T) where {T<:AbstractFloat} = sqrt(x^2 + eps_hyper^2)

"""
    compute_cz0hr_mapping!(Ri_reg, D_alpha, C_scale, Ri_raw, d2Ri_dz2, dz, config)

In-place, zero-allocation C-Z0HR mapping operator designed for inner-loop model integration.
Mutates pre-allocated vectors `Ri_reg`, `D_alpha`, and `C_scale`.
"""
function compute_cz0hr_mapping!(
    Ri_reg::AbstractVector{T},
    D_alpha::AbstractVector{T},
    C_scale::AbstractVector{T},
    Ri_raw::AbstractVector{T},
    d2Ri_dz2::AbstractVector{T},
    dz::T,
    config::CZ0HRConfig{T}
) where {T<:AbstractFloat}
    @inbounds for i in eachindex(Ri_raw)
        x = Ri_raw[i] - config.Ri_c
        Phi = Phi_eps(x, config.eps_hyper)
        
        # Discretized Spatial Curvature Scale C
        C = T(0.5) * abs(d2Ri_dz2[i]) * (dz^2) + config.epsilon_c * config.Ri_c
        C_scale[i] = C

        # Regularized Richardson Number Mapping
        num_reg = Phi + C
        den_reg = Phi + (one(T) + config.alpha) * C
        Ri_reg[i] = config.Ri_c + x * (num_reg / den_reg)

        # Partial Slope Evaluation D_α^{(R)} = ∂Ri_g^reg / ∂Ri_g |_C
        num_D = Phi^2 + T(2.0) * (one(T) + config.alpha) * C * Phi + (one(T) + config.alpha) * (C^2)
        D_alpha[i] = num_D / (den_reg^2)
    end
    return nothing
end

"""
    plot_cz0hr_diagnostics(z, Ri_raw, Ri_reg, Km_raw, Km_reg_Cinf)

Visualizes singularity regularization and C^∞ residual transport tail using Plots.jl.
"""
function plot_cz0hr_diagnostics(z, Ri_raw, Ri_reg, Km_raw, Km_reg_Cinf)
    if !@isdefined Plots
        return nothing
    end
    gr()
    p1 = plot(Ri_raw, z, xscale=:log10, label="Ri_raw", xlabel="Richardson Number (log)", ylabel="z (m)", title="Singularity Regularization")
    plot!(p1, Ri_reg, z, label="Ri_reg", lw=2)
    
    p2 = plot(Km_raw, z, label="Km (Raw C⁰)", xlabel="Km (m²/s)", ylabel="z (m)", title="Eddy Diffusivity Recovery")
    plot!(p2, Km_reg_Cinf, z, label="Km (C^∞ Soft-Plus)", lw=2)
    
    full_plot = plot(p1, p2, layout=(1, 2), size=(900, 450), margin=5Plots.mm)
    savefig(full_plot, "/workspace/scratch/cz0hr_llj_benchmark_v3.png")
    return full_plot
end

"""
    benchmark_cz0hr_llj_profile_v3()

Driver function for benchmarking the zero-allocation mutating operator in a 1D Low-Level Jet.
"""
function benchmark_cz0hr_llj_profile_v3()
    config = CZ0HRConfig()
    
    Nz = 300
    z = collect(range(1.0, 150.0, length=Nz))
    dz = z[2] - z[1]

    # Low-Level Jet Kinematics
    z_J = 50.0   # Jet core height [m]
    U_J = 10.0   # Jet core speed [m/s]

    U = [U_J * (zi / z_J) * exp(1.0 - zi / z_J) for zi in z]
    dU_dz = [(U_J / z_J) * (1.0 - zi / z_J) * exp(1.0 - zi / z_J) for zi in z]
    S_eff = [sqrt(du^2 + 1e-6) for du in dU_dz]

    g = 9.81
    theta_0 = 290.0
    dtheta_dz = [0.015 + 0.01 * exp(-zi / 30.0) for zi in z]
    N2 = [(g / theta_0) * dt for dt in dtheta_dz]

    Ri_raw = [N2[i] / (S_eff[i]^2) for i in 1:Nz]

    d2Ri_dz2 = zeros(Float64, Nz)
    for i in 2:(Nz-1)
        d2Ri_dz2[i] = (Ri_raw[i+1] - 2.0 * Ri_raw[i] + Ri_raw[i-1]) / (dz^2)
    end
    d2Ri_dz2[1] = d2Ri_dz2[2]
    d2Ri_dz2[end] = d2Ri_dz2[end-1]

    # Pre-allocated arrays for zero-allocation mutating call
    Ri_reg = zeros(Float64, Nz)
    D_alpha = zeros(Float64, Nz)
    C_scale = zeros(Float64, Nz)

    # Execute in-place zero-allocation operator
    compute_cz0hr_mapping!(Ri_reg, D_alpha, C_scale, Ri_raw, d2Ri_dz2, dz, config)

    # Exchange coefficients
    l_m = [15.0 * (zi / (zi + 20.0)) for zi in z]
    
    # Raw C^0 stability function
    Sm_raw_C0 = [max(0.0, (1.0 - Ri_raw[i] / config.Ri_c))^2 for i in 1:Nz]
    
    # End-to-end C^∞ Soft-Plus stability function
    g_raw_reg = 1.0 .- (Ri_reg ./ config.Ri_c)
    Sm_reg_Cinf = [(0.5 * (g_i + sqrt(g_i^2 + config.eps_hyper^2)))^2 for g_i in g_raw_reg]

    Km_raw = [(l_m[i]^2) * S_eff[i] * Sm_raw_C0[i] for i in 1:Nz]
    Km_reg_Cinf = [(l_m[i]^2) * S_eff[i] * Sm_reg_Cinf[i] for i in 1:Nz]

    # Render diagnostics
    plot_cz0hr_diagnostics(z, Ri_raw, Ri_reg, Km_raw, Km_reg_Cinf)

    return z, U, S_eff, Ri_raw, Ri_reg, D_alpha, Km_raw, Km_reg_Cinf
end

benchmark_cz0hr_llj_profile_v3()
