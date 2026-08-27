#!/usr/bin/env julia
# src/gspt_saddle_node_trajectory.jl

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
    fast_manifold_all(e, S, p)

Single-pass evaluator returning `(F, Fe, FS, Fee, FeS)`.
Note: The analytical derivatives assume S > 1e-4 (where dS_eff/dS = 1).
"""
function fast_manifold_all(e::Float64, S::Float64, p::SBLParams)
    S_eff = max(S, 1e-4)
    denom_S = S_eff^2 + p.beta * p.N2
    D = denom_S / (S_eff^2)
    denom_e = e^2 + p.delta_reg^2

    B0 = p.B0_max * (e^2 / denom_e)
    F = p.l0 * e * S_eff^2 - B0 - (e^3) / (p.l0 * D)

    B_e = 2.0 * p.B0_max * e * p.delta_reg^2 / (denom_e^2)
    B_ee = 2.0 * p.B0_max * p.delta_reg^2 * (p.delta_reg^2 - 3.0 * e^2) / (denom_e^3)

    Fe = p.l0 * S_eff^2 - B_e - (3.0 * e^2) / (p.l0 * D)
    FS = 2.0 * p.l0 * e * S_eff - (2.0 * p.beta * p.N2 * e^3 * S_eff) / (p.l0 * (denom_S^2))
    Fee = -B_ee - (6.0 * e) / (p.l0 * D)
    FeS = 2.0 * p.l0 * S_eff - (6.0 * p.beta * p.N2 * e^2 * S_eff) / (p.l0 * (denom_S^2))

    return F, Fe, FS, Fee, FeS
end

"""
    solve_asymptotic_fold(p)

Computes classical fold boundary as δ -> 0, with fallback for N2=0 or β=0.
"""
function solve_asymptotic_fold(p::SBLParams)
    target = 27.0 * p.B0_max^2 / (4.0 * p.l0^4)
    a = p.beta * p.N2

    # Degenerate unstratified/neutral case (a = 0)
    if a < 1e-12
        S_fold = target^(1/6)
        e_fold = 3.0 * p.B0_max / (2.0 * p.l0 * S_fold^2)
        return S_fold, e_fold
    end

    p_c = -a^2 / 3.0
    q_c = 2.0 * a^3 / 27.0 - target

    arg = (3.0 * q_c) / (2.0 * p_c) * sqrt(-3.0 / p_c)
    if abs(arg) > 1.0
        @warn "Cardano discriminant out of range (|arg| = $(abs(arg)) > 1)."
    end

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
    solve_regularized_fold(p; tol=1e-10, max_iter=50)

2D Newton-Raphson solver initialized dynamically from the asymptotic fold.
"""
function solve_regularized_fold(p::SBLParams; tol::Float64=1e-10, max_iter::Int=50)
    S_asym, e_asym = solve_asymptotic_fold(p)
    x = [e_asym, S_asym] # Dynamic initialization

    for iter in 1:max_iter
        e, S = x[1], x[2]
        F, Fe, FS, Fee, FeS = fast_manifold_all(e, S, p)
        res = [F, Fe]

        if norm(res, Inf) < tol
            if x[1] <= 0 || x[2] <= 0
                @warn "Newton solver converged to non-physical domain: e=$(x[1]), S=$(x[2])"
                return NaN, NaN
            end
            return S, e
        end

        J = [Fe FS; Fee FeS]
        if abs(det(J)) < 1e-14
            @warn "Singular Jacobian encounter in Newton fold solver."
            return NaN, NaN
        end
        x -= J \ res
    end

    @warn "Newton solver failed to converge within $max_iter iterations."
    return NaN, NaN
end

"""
    simulate_sbl_dynamics(; t_max=80.0, dt=0.001)
"""
function simulate_sbl_dynamics(; t_max::Float64=80.0, dt::Float64=0.001)
    p = SBLParams()
    t_arr = collect(0.0:dt:t_max)
    n_steps = length(t_arr)

    e_arr, S_arr = zeros(n_steps), zeros(n_steps)
    lambda_fast, norm_dist = zeros(n_steps), zeros(n_steps)

    e_arr[1], S_arr[1] = 0.5, 1.0

    # Diagnostics tracking
    floor_episodes = 0
    in_floor = false
    max_euler_ratio = 0.0

    for i in 2:n_steps
        e, S = e_arr[i-1], S_arr[i-1]
        F, Fe, FS, _, _ = fast_manifold_all(e, S, p)

        # Explicit Euler stability CFL check: dt < 2 * epsilon / |Fe|
        current_ratio = dt * abs(Fe) / (2.0 * p.epsilon)
        max_euler_ratio = max(max_euler_ratio, current_ratio)

        de_dt = (1.0 / p.epsilon) * F
        dS_dt = p.G0 - p.gamma_s * e * S - p.r_s * S

        e_next = e + de_dt * dt
        S_next = S + dS_dt * dt

        if e_next <= p.e_floor
            e_next = p.e_floor
            if !in_floor
                floor_episodes += 1
                in_floor = true
            end
        else
            in_floor = false
        end

        e_arr[i] = e_next
        S_arr[i] = max(S_next, 0.0)

        # Post-step fast subsystem diagnostics
        F_i, Fe_i, FS_i, _, _ = fast_manifold_all(e_arr[i], S_arr[i], p)
        lambda_fast[i] = (1.0 / p.epsilon) * Fe_i
        norm_dist[i] = abs(F_i) / max(1.0, abs(Fe_i) * e_arr[i] + abs(FS_i) * S_arr[i])
    end

    S_fold_asym, e_fold_asym = solve_asymptotic_fold(p)
    S_fold_reg, e_fold_reg = solve_regularized_fold(p)

    # Consistency Assertion: Regularized fold must match asymptotic fold as δ -> 0
    @assert abs(S_fold_reg - S_fold_asym) < 1e-2 "Fold solver consistency check failed."

    if max_euler_ratio > 1.0
        @warn @sprintf("Forward Euler step size exceeds stability limit (Max dt*|λ_fast|/2 = %.2f > 1.0).", max_euler_ratio)
    end

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
    @printf("Positivity Floor Episodes   : %d discrete event(s)\n", floor_episodes)
    @printf("Max Euler Stability Ratio  : %.4f (Limit < 1.0)\n", max_euler_ratio)
    println("="^80)

    println("\nTrajectory Diagnostic Log (GSPT Verification):")
    println("-"^80)
    @printf("%-8s | %-10s | %-10s | %-12s | %-12s | %-15s\n",
        "Time (s)", "Shear S", "TKE e", "λ_fast", "d_C*", "Hyperbolicity State")
    println("-"^80)
    step_stride = max(1, round(Int, n_steps / 15))
    for i in 1:step_stride:n_steps
        # Fast timescale threshold: |λ_fast| > 0.1 / ε (dimensionless timescale separation)
        state = lambda_fast[i] < -0.1 / p.epsilon ? "Attracting" :
                (lambda_fast[i] > 0.1 / p.epsilon ? "Repelling" : "Fold Trans.")
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