#!/usr/bin/env julia
# src/gspt_saddle_node_trajectory.jl
# version 0.2.0
using LinearAlgebra
using Printf

Base.@kwdef struct SBLParams
    epsilon::Float64 = 0.05    # Fast time-scale parameter (TKE relaxation)
    l0::Float64 = 1.0     # Mixing length scale
    beta::Float64 = 5.0     # Stability coefficient
    N2::Float64 = 0.1     # Stratification (Brunt-Vaisala frequency squared)
    B0_max::Float64 = 0.05    # Base buoyancy flux
    delta_reg::Float64 = 0.01    # Buoyancy flux regularization threshold
    G0::Float64 = 0.3     # Geostrophic wind shear forcing
    gamma_s::Float64 = 1.8     # Turbulent shear destruction gain
    r_s::Float64 = 0.15    # Background linear shear relaxation
    e_floor::Float64 = 1e-4    # Positivity protection threshold
end

"""
    fast_manifold_F(e, S, p)

Evaluates the fast subsystem reduced vector field residual F(e, S; δ) where:
    f_fast = (1/ε) * F(e, S; δ)
"""
function fast_manifold_F(e::Float64, S::Float64, p::SBLParams)
    S_eff = max(S, 1e-4)
    D = 1.0 + p.beta * p.N2 / (S_eff^2)
    B0 = p.B0_max * (e^2 / (e^2 + p.delta_reg^2))
    return p.l0 * e * S_eff^2 - B0 - (e^3) / (p.l0 * D)
end

"""
    fast_manifold_jac(e, S, p)

Computes analytical partial derivatives F_e, F_S, F_ee, and F_eS for stability and fold solvers.
"""
function fast_manifold_jac(e::Float64, S::Float64, p::SBLParams)
    S_eff = max(S, 1e-4)
    D = 1.0 + p.beta * p.N2 / (S_eff^2)
    denom_e = e^2 + p.delta_reg^2
    denom_S = S_eff^2 + p.beta * p.N2

    B_e = 2.0 * p.B0_max * e * p.delta_reg^2 / (denom_e^2)
    B_ee = 2.0 * p.B0_max * p.delta_reg^2 * (p.delta_reg^2 - 3.0 * e^2) / (denom_e^3)

    Fe = p.l0 * S_eff^2 - B_e - (3.0 * e^2) / (p.l0 * D)
    FS = 2.0 * p.l0 * e * S_eff - (2.0 * p.beta * p.N2 * e^3 * S_eff) / (p.l0 * (denom_S^2))
    Fee = -B_ee - (6.0 * e) / (p.l0 * D)
    FeS = 2.0 * p.l0 * S_eff - (6.0 * p.beta * p.N2 * e^2 * S_eff) / (p.l0 * (denom_S^2))

    return Fe, FS, Fee, FeS
end

"""
    solve_regularized_fold(p; tol=1e-10, max_iter=50)

Exact saddle-node fold solver for the regularized system (δ > 0).
Solves the simultaneous nonlinear algebraic system:
    F(e, S; δ) = 0
    F_e(e, S; δ) = 0
using a 2D Newton-Raphson scheme.
"""
function solve_regularized_fold(p::SBLParams; tol::Float64=1e-10, max_iter::Int=50)
    x = [0.47, 0.40] # Initial guess (e, S) near the fold region
    for _ in 1:max_iter
        e, S = x[1], x[2]
        F = fast_manifold_F(e, S, p)
        Fe, FS, Fee, FeS = fast_manifold_jac(e, S, p)

        res = [F, Fe]
        if norm(res, Inf) < tol
            return S, e
        end

        J = [Fe FS; Fee FeS]
        x -= J \ res
    end
    return NaN, NaN
end

"""
    solve_asymptotic_fold(p)

Computes the classical cubic saddle-node fold boundary as δ -> 0.
"""
function solve_asymptotic_fold(p::SBLParams)
    target = 27.0 * p.B0_max^2 / (4.0 * p.l0^4)
    a = p.beta * p.N2
    p_c = -a^2 / 3.0
    q_c = 2.0 * a^3 / 27.0 - target

    arg = (3.0 * q_c) / (2.0 * p_c) * sqrt(-3.0 / p_c)
    phi = acos(clamp(arg, -1.0, 1.0))
    R = 2.0 * sqrt(-p_c / 3.0)

    x_max = max(
        R * cos(phi / 3.0) - a / 3.0,
        R * cos((phi + 2.0 * pi) / 3.0) - a / 3.0,
        R * cos((phi + 4.0 * pi) / 3.0) - a / 3.0
    )
    S_fold = sqrt(x_max)
    e_fold = 3.0 * p.B0_max / (2.0 * p.l0 * S_fold^2)
    return S_fold, e_fold
end

"""
    simulate_sbl_dynamics(; t_max=80.0, dt=0.001)

Executes slow-fast trajectory integration alongside fast eigenvalue analysis,
manifold distance diagnostics, and dynamic fold departure tracking.
"""
function simulate_sbl_dynamics(; t_max::Float64=80.0, dt::Float64=0.001)
    p = SBLParams()
    t_arr = collect(0.0:dt:t_max)
    n_steps = length(t_arr)

    e_arr = zeros(n_steps)
    S_arr = zeros(n_steps)
    lambda_fast = zeros(n_steps)
    norm_dist = zeros(n_steps)

    # Initial conditions on the attracting branch
    e_arr[1], S_arr[1] = 0.5, 1.0
    floor_hits = 0

    for i in 2:n_steps
        e, S = e_arr[i-1], S_arr[i-1]

        # Fast variable update
        F = fast_manifold_F(e, S, p)
        de_dt = (1.0 / p.epsilon) * F

        # Slow variable update
        dS_dt = p.G0 - p.gamma_s * e * S - p.r_s * S

        e_next = e + de_dt * dt
        S_next = S + dS_dt * dt

        # Track positivity floor activations
        if e_next < p.e_floor
            e_next = p.e_floor
            floor_hits += 1
        end

        e_arr[i] = e_next
        S_arr[i] = max(S_next, 0.0)

        # Fast subsystem diagnostics
        Fe, FS, _, _ = fast_manifold_jac(e_arr[i], S_arr[i], p)
        lambda_fast[i] = (1.0 / p.epsilon) * Fe
        F_i = fast_manifold_F(e_arr[i], S_arr[i], p)

        # Normalized distance d_C* = |F| / max(1, |F_e|e + |F_S|S)
        norm_dist[i] = abs(F_i) / max(1.0, abs(Fe) * e_arr[i] + abs(FS) * S_arr[i])
    end

    # Static fold calculations
    S_fold_asym, e_fold_asym = solve_asymptotic_fold(p)
    S_fold_reg, e_fold_reg = solve_regularized_fold(p)

    # Identify dynamic loss of normal hyperbolicity (first λ_fast > 0 transition)
    dep_idx = findfirst(i -> lambda_fast[i] > 0.0, 1:n_steps)
    S_dep = dep_idx !== nothing ? S_arr[dep_idx] : NaN
    e_dep = dep_idx !== nothing ? e_arr[dep_idx] : NaN

    println("="^80)
    println("                      GSPT MANIFOLD & BIFURCATION AUDIT RESULTS")
    println("="^80)
    @printf("Asymptotic Fold (δ = 0.00)  : S_fold = %8.5f s^-1, e_fold = %8.5f m^2 s^-2\n", S_fold_asym, e_fold_asym)
    @printf("Regularized Fold (δ = 0.01) : S_fold = %8.5f s^-1, e_fold = %8.5f m^2 s^-2\n", S_fold_reg, e_fold_reg)
    if dep_idx !== nothing
        @printf("Dynamic Departure (λ_fast>0): S_dep  = %8.5f s^-1, e_dep  = %8.5f m^2 s^-2 (t = %.2f s)\n",
            S_dep, e_dep, t_arr[dep_idx])
    end
    @printf("Numerical Floor Activations : %d / %d integration steps (%.2f%%)\n",
        floor_hits, n_steps, (floor_hits / n_steps) * 100)
    println("="^80)

    # Trajectory Diagnostics Sample Log
    println("\nTrajectory Diagnostic Log (GSPT Verification):")
    println("-"^80)
    @printf("%-8s | %-10s | %-10s | %-12s | %-12s | %-15s\n",
        "Time (s)", "Shear S", "TKE e", "λ_fast", "d_C*", "Hyperbolicity State")
    println("-"^80)
    step_stride = max(1, round(Int, n_steps / 15))
    for i in 1:step_stride:n_steps
        state = lambda_fast[i] < -0.1 ? "Attracting" : (lambda_fast[i] > 0.1 ? "Repelling" : "Fold Trans.")
        @printf("%8.2f | %10.4f | %10.4f | %12.4e | %12.4e | %-15s\n",
            t_arr[i], S_arr[i], e_arr[i], lambda_fast[i], norm_dist[i], state)
    end
    println("="^80)

    return S_fold_reg, e_fold_reg, S_arr, e_arr, lambda_fast, norm_dist
end

function main()
    simulate_sbl_dynamics()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end