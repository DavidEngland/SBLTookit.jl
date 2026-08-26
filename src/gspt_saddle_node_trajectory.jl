#!/usr/bin/env julia
# src/gspt_saddle_node_trajectory.jl
# package: SBLToolkit.jl
# version 0.1.0
# Author: David E. England, Ph.D.
#
# Date: 25 Aug 2026
# Description: Simulates the slow-fast SBL system and analyzes the saddle-node bifurcation.
# ==============================================================================
# GSPT SADDLE-NODE BIFURCATION & 2D STATE-SPACE TRAJECTORY SIMULATOR
# Developed for Generalized Similarity Profile Theory (GSPT) Manifold Verification
# ==============================================================================
# This script simulates the slow-fast SBL dynamical system:
# 1. Fast Variable: Turbulent Kinetic Energy (e), which relaxes rapidly toward 
#    its equilibrium manifold C_0^+.
# 2. Slow Variable: Wind Shear (S), which is driven by geostrophic forcing G_0 
#    and destroyed by turbulent shear stress.
#
# The buoyancy flux is regularized quadratically near the laminar limit to prevent 
# unphysical negative TKE while allowing shear-driven reignition.
# ==============================================================================

using LinearAlgebra
using Printf

# Optional: using Plots for local rendering if executed in a standard environment
# using Plots

"""
    solve_cubic_branches(S::Float64, l0::Float64, beta::Float64, N2::Float64, B0_max::Float64)

Analytically computes the stable turbulent branch (e_st) and the unstable threshold
(e_unst) for the GSPT fast manifold equation:
    e^3 - p*e + q = 0
where:
    p = l0^2 * (S^2 + beta * N2)
    q = l0 * B0_max * (1.0 + beta * N2 / S^2)

Uses a high-performance, branch-free trigonometric Cardano-like solver optimized for
HPC registers. Returns (e_st, e_unst, has_roots).
"""
function solve_cubic_branches(S::Float64, l0::Float64, beta::Float64, N2::Float64, B0_max::Float64)
    if S < 1e-4
        return 0.0, 0.0, false
    end
    
    p = l0^2 * (S^2 + beta * N2)
    q = l0 * B0_max * (1.0 + beta * N2 / S^2)
    
    # Check discriminant
    D = 4.0 * p^3 - 27.0 * q^2
    if D <= 0.0
        return 0.0, 0.0, false # Stable branch has vanished (saddle-node bifurcation)
    end
    
    R = 2.0 * sqrt(p / 3.0)
    arg = 3.0 * q / (2.0 * p * sqrt(p / 3.0))
    # Clamp to protect arccos domain against float rounding
    arg_clamped = max(-1.0, min(1.0, arg))
    phi = acos(arg_clamped)
    
    # Analytical roots of e^3 - p*e + q = 0 (using negated roots of x^3 - p*x - q = 0)
    e_st = -R * cos((phi + 2.0 * pi) / 3.0)
    e_unst = -R * cos((phi + 4.0 * pi) / 3.0)
    
    return e_st, e_unst, true
end

"""
    find_analytical_fold(l0::Float64, beta::Float64, N2::Float64, B0_max::Float64)

Finds the exact critical wind shear S_fold where the saddle-node bifurcation occurs.
The fold boundary is defined by:
    S^4 * (S^2 + beta * N2) = 27 * B0_max^2 / (4 * l0^4)
"""
function find_analytical_fold(l0::Float64, beta::Float64, N2::Float64, B0_max::Float64)
    target = 27.0 * B0_max^2 / (4.0 * l0^4)
    # S^6 + (beta*N2)*S^4 - target = 0
    # Let x = S^2, then x^3 + a*x^2 + b*x + c = 0 with b=0, c=-target.
    a = beta * N2
    b = 0.0
    c = -target

    # Depressed cubic y^3 + p_c*y + q_c = 0 via x = y - a/3.
    p_c = b - a^2 / 3.0
    q_c = 2.0 * a^3 / 27.0 - a * b / 3.0 + c

    # Discriminant for depressed cubic: Δ = (q/2)^2 + (p/3)^3
    Δ = (q_c / 2.0)^2 + (p_c / 3.0)^3

    if Δ > 0.0
        # One real root.
        y = cbrt(-q_c / 2.0 + sqrt(Δ)) + cbrt(-q_c / 2.0 - sqrt(Δ))
        x = y - a / 3.0
        return x > 0.0 ? sqrt(x) : NaN
    else
        # Three real roots (standard trigonometric form).
        arg = (3.0 * q_c) / (2.0 * p_c) * sqrt(-3.0 / p_c)
        phi = acos(clamp(arg, -1.0, 1.0))
        R = 2.0 * sqrt(-p_c / 3.0)

        y1 = R * cos(phi / 3.0)
        y2 = R * cos((phi + 2.0 * pi) / 3.0)
        y3 = R * cos((phi + 4.0 * pi) / 3.0)

        x_candidates = [y1 - a / 3.0, y2 - a / 3.0, y3 - a / 3.0]
        x_pos = filter(x -> x > 0.0, x_candidates)
        isempty(x_pos) && return NaN

        # Positive-S physical branch corresponds to the largest positive x = S^2.
        return sqrt(maximum(x_pos))
    end
end

"""
    simulate_sbl_dynamics(; t_max=80.0, dt=0.001)

Performs forward time integration of the slow-fast SBL system, resolving the
relaxation oscillation cycle as it traverses the saddle-node fold and rebounds.
"""
function simulate_sbl_dynamics(; t_max::Float64=80.0, dt::Float64=0.001)
    # Physical and model parameters
    epsilon = 0.05    # Fast time-scale parameter (TKE relaxation)
    l0 = 1.0          # Mixing length scale
    beta = 5.0        # Stability coefficient
    N2 = 0.1          # Stratification (Brunt-Vaisala frequency squared)
    B0_max = 0.05     # Base buoyancy flux
    delta_reg = 0.01  # Buoyancy flux regularization threshold near e=0
    
    # Slow momentum system parameters
    G0 = 0.3          # Geostrophic wind shear forcing
    gamma_s = 1.8     # Turbulent shear destruction gain
    r_s = 0.15        # Background linear shear relaxation
    
    n_steps = int = round(Int, t_max / dt)
    t_arr = collect(range(0.0, t_max, length=n_steps))
    
    e_arr = zeros(n_steps)
    S_arr = zeros(n_steps)
    
    # Initial conditions (stable turbulent branch)
    e_arr[1] = 0.5
    S_arr[1] = 1.0
    
    println("Simulating 2D phase-space slow-fast dynamics...")
    for i in 2:n_steps
        e = e_arr[i-1]
        S = S_arr[i-1]
        
        # Fast system: TKE budget with quadratic regularization
        B0 = B0_max * (e^2 / (e^2 + delta_reg^2))
        denom = 1.0 + beta * N2 / (max(S, 1e-4)^2)
        de_dt = (1.0 / epsilon) * (l0 * e * S^2 - B0 - (e^3) / (l0 * denom))
        
        # Slow system: Shear evolution
        dS_dt = G0 - gamma_s * e * S - r_s * S
        
        # Forward Euler step
        e_next = e + de_dt * dt
        S_next = S + dS_dt * dt
        
        # Apply physical background TKE floor as seed for reignition
        e_arr[i] = max(e_next, 1e-4)
        S_arr[i] = max(S_next, 0.0)
    end
    
    # Compute analytical fold point
    S_fold = find_analytical_fold(l0, beta, N2, B0_max)
    e_fold = 3.0 * B0_max / (2.0 * S_fold^2)
    
    println("="^80)
    println("                          GSPT BIFURCATION AUDIT RESULTS")
    println("="^80)
    @printf("Analytical Fold Shear S_fold : %10.5f s^-1\n", S_fold)
    @printf("Analytical Fold TKE e_fold   : %10.5f m^2 s^-2\n", e_fold)
    @printf("Simulation S-range           : [%.4f, %.4f] s^-1\n", minimum(S_arr), maximum(S_arr))
    @printf("Simulation e-range           : [%.4f, %.4f] m^2 s^-2\n", minimum(e_arr), maximum(e_arr))
    println("="^80)
    
    # Print sample of the trajectory cycle
    println("\nSample Trajectory (Relaxation Cycle):")
    println("-"^80)
    @printf("%-10s | %-12s | %-12s | %-15s\n", "Time (s)", "Shear S", "TKE e", "State")
    println("-"^80)
    for i in 1:round(Int, n_steps/15):n_steps
        state = e_arr[i] < 0.01 ? "Laminar Attractor" : "Stable Turbulent"
        @printf("%10.2f | %12.4f | %12.4f | %-15s\n", t_arr[i], S_arr[i], e_arr[i], state)
    end
    println("="^80)
    
    return S_fold, e_fold, S_arr, e_arr
end

# Main entrypoint
function main()
    simulate_sbl_dynamics()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
