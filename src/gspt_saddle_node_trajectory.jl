#!/usr/bin/env julia
# src/gspt_saddle_node_trajectory.jl

using LinearAlgebra
using Printf

Base.@kwdef struct SBLParams
    epsilon::Float64 = 0.05    # Fast time-scale parameter
    l0::Float64 = 1.0     # Mixing length scale
    beta::Float64 = 5.0     # Stability coefficient
    N2::Float64 = 0.1     # Stratification (Brunt-Vaisala frequency squared)
    B0_max::Float64 = 0.05    # Base buoyancy flux
    delta_reg::Float64 = 0.01    # Buoyancy flux regularization threshold
    eta_S::Float64 = 1e-4    # C^infinity shear smoothing parameter
    G0::Float64 = 0.3     # Geostrophic wind shear forcing
    gamma_s::Float64 = 1.8     # Turbulent shear destruction gain
    r_s::Float64 = 0.15    # Background linear shear relaxation
    e_floor::Float64 = 1e-4    # Positivity protection threshold
end

"""
    SBLParams(p::SBLParams; kwargs...)

Constructs a copy of `p` with overridden parameter values.
"""
function SBLParams(p::SBLParams; kwargs...)
    kw = Dict(kwargs...)
    return SBLParams(
        get(kw, :epsilon, p.epsilon),
        get(kw, :l0, p.l0),
        get(kw, :beta, p.beta),
        get(kw, :N2, p.N2),
        get(kw, :B0_max, p.B0_max),
        get(kw, :delta_reg, p.delta_reg),
        get(kw, :eta_S, p.eta_S),
        get(kw, :G0, p.G0),
        get(kw, :gamma_s, p.gamma_s),
        get(kw, :r_s, p.r_s),
        get(kw, :e_floor, p.e_floor)
    )
end

"""
    fast_manifold_all(e, S, p)

Single-pass C^infinity evaluator returning `(F, Fe, FS, Fee, FeS)`.
Uses smooth regularization S_eff = sqrt(S^2 + eta_S^2).
"""
function fast_manifold_all(e::Float64, S::Float64, p::SBLParams)
    S_eff = sqrt(S^2 + p.eta_S^2)
    dS_eff_dS = S / S_eff

    a = p.beta * p.N2
    denom_S = S_eff^2 + a
    D = denom_S / (S_eff^2)
    denom_e = e^2 + p.delta_reg^2

    B0 = p.B0_max * (e^2 / denom_e)
    F = p.l0 * e * S_eff^2 - B0 - (e^3) / (p.l0 * D)

    B_e = 2.0 * p.B0_max * e * p.delta_reg^2 / (denom_e^2)
    B_ee = 2.0 * p.B0_max * p.delta_reg^2 * (p.delta_reg^2 - 3.0 * e^2) / (denom_e^3)

    Fe = p.l0 * S_eff^2 - B_e - (3.0 * e^2) / (p.l0 * D)

    dF_dS_eff = 2.0 * p.l0 * e * S_eff - (2.0 * a * e^3 * S_eff) / (p.l0 * (denom_S^2))
    FS = dF_dS_eff * dS_eff_dS

    Fee = -B_ee - (6.0 * e) / (p.l0 * D)

    dFe_dS_eff = 2.0 * p.l0 * S_eff - (6.0 * a * e^2 * S_eff) / (p.l0 * (denom_S^2))
    FeS = dFe_dS_eff * dS_eff_dS

    return F, Fe, FS, Fee, FeS
end

"""
    solve_asymptotic_fold(p; eta_S=p.eta_S, tol=1e-12, max_iter=30)

Solves the asymptotic fold (δ = 0). Computes x = S_eff^2 via x^3 + a*x^2 - C = 0,
then recovers S_fold = sqrt(x - eta_S^2). Initialized at x0 = cbrt(C).
"""
function solve_asymptotic_fold(p::SBLParams; eta_S::Float64=p.eta_S, tol::Float64=1e-12, max_iter::Int=30)
    @assert p.beta >= 0.0 && p.N2 >= 0.0 "Physical parameters beta and N2 must be non-negative."
    a = p.beta * p.N2
    C = 27.0 * p.B0_max^2 / (4.0 * p.l0^4)

    x = cbrt(C) # Guaranteed upper bound initialization: f(cbrt(C)) = a*C^(2/3) >= 0
    for iter in 1:max_iter
        fx = x^3 + a * x^2 - C
        fpx = 3.0 * x^2 + 2.0 * a * x
        if abs(fpx) < 1e-14
            @warn "Near-zero derivative encounter in 1D fold solver at iteration $iter."
            return NaN, NaN
        end
        dx = fx / fpx
        x -= dx
        if abs(dx) < tol
            S_eff2 = max(x, 0.0)
            S_fold2 = S_eff2 - eta_S^2
            if S_fold2 < 0.0
                @warn "Square of S_fold is negative under eta_S correction."
                return NaN, NaN
            end
            S_fold = sqrt(S_fold2)
            e_fold = 3.0 * p.B0_max / (2.0 * p.l0 * S_eff2)
            return S_fold, e_fold
        end
    end

    @warn "1D Asymptotic fold solver failed to converge within $max_iter iterations."
    return NaN, NaN
end

"""
    solve_regularized_fold(p; tol=1e-10, max_iter=50)

2D Newton-Raphson solver for the regularized fold (F = 0, Fe = 0).
"""
function solve_regularized_fold(p::SBLParams; tol::Float64=1e-10, max_iter::Int=50)
    S_asym, e_asym = solve_asymptotic_fold(p)
    if isnan(S_asym)
        return NaN, NaN
    end
    x = [e_asym, S_asym]

    for iter in 1:max_iter
        e, S = x[1], x[2]
        F, Fe, FS, Fee, FeS = fast_manifold_all(e, S, p)
        res = [F, Fe]

        if norm(res, Inf) < tol
            if x[1] <= 0 || x[2] <= 0
                @warn "Newton solver converged outside physical domain (iter $iter): e=$(x[1]), S=$(x[2])"
                return NaN, NaN
            end
            return S, e
        end

        J = [Fe FS; Fee FeS]
        if abs(det(J)) < 1e-14
            @warn "Singular Jacobian encounter in Newton fold solver at iteration $iter."
            return NaN, NaN
        end
        x -= J \ res
    end

    @warn "Newton fold solver failed to converge within $max_iter iterations."
    return NaN, NaN
end

"""
    initialize_on_manifold(S0, p; grid_pts=200, e_max=5.0)

Scans e-space to select the root of F(e, S0) = 0 on the attracting sheet C0+ (Fe < 0).
"""
function initialize_on_manifold(S0::Float64, p::SBLParams; grid_pts::Int=200, e_max::Float64=5.0)
    e_grid = range(1e-3, e_max, length=grid_pts)
    candidates = Float64[]

    for e_guess in e_grid
        e = e_guess
        for _ in 1:20
            F, Fe, _, _, _ = fast_manifold_all(e, S0, p)
            if abs(F) < 1e-11
                if Fe < 0.0 && !(any(c -> abs(c - e) < 1e-5, candidates))
                    push!(candidates, e)
                end
                break
            end
            if abs(Fe) < 1e-14
                ;
                break;
            end
            e -= F / Fe
        end
    end

    if isempty(candidates)
        @warn "No attracting branch root (Fe < 0) found for S0 = $S0. Defaulting to e = 0.8."
        return 0.8
    end

    return maximum(candidates)
end

"""
    audit_fold_nondegeneracy(S_fold, e_fold, p)

Verifies saddle-node transversality (FS != 0) and nondegeneracy (Fee != 0).
"""
function audit_fold_nondegeneracy(S_fold::Float64, e_fold::Float64, p::SBLParams)
    F, Fe, FS, Fee, FeS = fast_manifold_all(e_fold, S_fold, p)
    detJ = -FS * Fee

    println("Saddle-Node Nondegeneracy & Transversality Audit:")
    @printf("  Fold Residual |F|         : %12.4e\n", abs(F))
    @printf("  Fold Residual |Fe|        : %12.4e\n", abs(Fe))
    @printf("  Transversality (FS)       : %12.5f (Requires FS != 0)\n", FS)
    @printf("  Nondegeneracy (Fee)       : %12.5f (Requires Fee != 0)\n", Fee)
    @printf("  Jacobian Det |det J_fold| : %12.5f\n", abs(detJ))
    println()
end

"""
    audit_2d_fold_convergence(p_base)

Evaluates joint regularization limit (δ, η_S) -> (0, 0) against the unregularized fold.
"""
function audit_2d_fold_convergence(p_base::SBLParams)
    S_classical, e_classical = solve_asymptotic_fold(p_base; eta_S=0.0)

    println("\n2D Parameter Regularization Audit ((δ, η_S) -> (0, 0)):")
    println("-"^75)
    @printf("%-10s | %-10s | %-12s | %-12s | %-12s | %-12s\n", "δ_reg", "η_S", "S_fold", "e_fold", "ΔS_fold", "Δe_fold")
    println("-"^75)

    test_pairs = [(1e-1, 1e-2), (5e-2, 5e-3), (1e-2, 1e-4), (1e-3, 1e-5), (1e-4, 1e-6)]
    last_dS = Inf
    monotone = true
    max_err = 0.0

    for (d, eta) in test_pairs
        p_test = SBLParams(p_base; delta_reg=d, eta_S=eta)
        S_reg, e_reg = solve_regularized_fold(p_test)
        dS = abs(S_reg - S_classical)
        de = abs(e_reg - e_classical)
        max_err = max(max_err, dS)

        if dS > last_dS + 1e-12
            monotone = false
        end
        last_dS = dS

        @printf("%10.1e | %10.1e | %12.5f | %12.5f | %12.4e | %12.4e\n", d, eta, S_reg, e_reg, dS, de)
    end
    println("-"^75)

    if max_err < 0.05
        println("VERIFIED: Fold converges to the unregularized reference limit within tolerance.")
    else
        println("VERIFICATION WARNING: Divergence detected during parameter sweep.")
    end
    println("Observed Monotonicity: ", monotone ? "Yes" : "No", "\n")
end

"""
    simulate_sbl_dynamics(; t_max=80.0, dt=0.001, init_on_manifold=true)
"""
function simulate_sbl_dynamics(; t_max::Float64=80.0, dt::Float64=0.001, init_on_manifold::Bool=true)
    p = SBLParams()
    t_arr = collect(0.0:dt:t_max)
    n_steps = length(t_arr)

    e_arr, S_arr = zeros(n_steps), zeros(n_steps)
    lambda_fast, manifold_residual = zeros(n_steps), zeros(n_steps)

    S_arr[1] = 1.0
    e_arr[1] = init_on_manifold ? initialize_on_manifold(S_arr[1], p) : 0.5

    floor_episodes = 0
    in_floor = false
    max_fast_euler_ratio = 0.0

    # GSPT Event Tracking Indices
    idx_fold_crossing = 0
    idx_repelling_onset = 0
    idx_dynamic_departure = 0

    S_fold_reg, e_fold_reg = solve_regularized_fold(p)

    for i in 2:n_steps
        e, S = e_arr[i-1], S_arr[i-1]
        F, Fe, FS, _, _ = fast_manifold_all(e, S, p)

        ratio = dt * (abs(Fe) / p.epsilon) / 2.0
        max_fast_euler_ratio = max(max_fast_euler_ratio, ratio)

        de_dt = (1.0 / p.epsilon) * F
        dS_dt = p.G0 - p.gamma_s * e * S - p.r_s * S

        e_next = e + de_dt * dt
        S_next = S + dS_dt * dt

        # Positivity safeguard
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

        F_i, Fe_i, FS_i, _, _ = fast_manifold_all(e_arr[i], S_arr[i], p)
        lambda_fast[i] = (1.0 / p.epsilon) * Fe_i
        res_i = abs(F_i) / max(1.0, abs(Fe_i) * e_arr[i] + abs(FS_i) * S_arr[i])
        manifold_residual[i] = res_i

        # Event 1: Geometric Fold Crossing (S <= S_fold)
        if idx_fold_crossing == 0 && S_arr[i] <= S_fold_reg
            idx_fold_crossing = i
        end

        # Event 2: Linear Repelling Onset (Fe > 0)
        if idx_repelling_onset == 0 && Fe_i > 0.0
            idx_repelling_onset = i
        end

        # Event 3: Dynamic Departure / Escape (residual threshold)
        if idx_repelling_onset > 0 && idx_dynamic_departure == 0 && res_i > 0.05
            idx_dynamic_departure = i
        end
    end

    S_fold_asym, e_fold_asym = solve_asymptotic_fold(p)

    if max_fast_euler_ratio > 1.0
        @warn @sprintf("Explicit Euler step size violates frozen-fast stability threshold (Max ratio = %.4f > 1.0).", max_fast_euler_ratio)
    end

    println("="^80)
    println("                      GSPT MANIFOLD & BIFURCATION AUDIT RESULTS")
    println("="^80)
    @printf("Asymptotic Fold (δ=0, η_S=%.0e) : S_fold = %8.5f s^-1, e_fold = %8.5f m^2 s^-2\n", p.eta_S, S_fold_asym, e_fold_asym)
    @printf("Regularized Fold (δ=%.2f, η_S=%.0e): S_fold = %8.5f s^-1, e_fold = %8.5f m^2 s^-2\n", p.delta_reg, p.eta_S, S_fold_reg, e_fold_reg)

    println("-"^80)
    println("GSPT Dynamic Trajectory Events:")
    if idx_fold_crossing > 0
        @printf("  1. Geometric Fold Crossing (S <= S_f) : t = %6.2f s | S = %.5f, e = %.5f\n",
            t_arr[idx_fold_crossing], S_arr[idx_fold_crossing], e_arr[idx_fold_crossing])
    else
        println("  1. Geometric Fold Crossing            : Not reached during integration.")
    end

    if idx_repelling_onset > 0
        @printf("  2. Linear Repelling Onset (λ > 0)     : t = %6.2f s | S = %.5f, e = %.5f\n",
            t_arr[idx_repelling_onset], S_arr[idx_repelling_onset], e_arr[idx_repelling_onset])
    else
        println("  2. Linear Repelling Onset             : Not triggered during integration.")
    end

    if idx_dynamic_departure > 0
        @printf("  3. Dynamic Trajectory Departure       : t = %6.2f s | S = %.5f, e = %.5f\n",
            t_arr[idx_dynamic_departure], S_arr[idx_dynamic_departure], e_arr[idx_dynamic_departure])
    else
        println("  3. Dynamic Trajectory Departure       : No boundary departure observed.")
    end

    println("-"^80)
    @printf("Positivity Floor Episodes               : %d discrete event(s)\n", floor_episodes)
    @printf("Max Fast Euler Ratio                    : %.4f (Frozen-fast limit < 1.0)\n", max_fast_euler_ratio)
    println("="^80)
    println()

    audit_fold_nondegeneracy(S_fold_reg, e_fold_reg, p)
    audit_2d_fold_convergence(p)

    println("Trajectory Diagnostic Log:")
    println("-"^80)
    @printf("%-8s | %-10s | %-12s | %-12s | %-12s | %-15s\n",
        "Time (s)", "Shear S", "TKE e", "λ_fast", "Residual", "Hyperbolicity State")
    println("-"^80)
    step_stride = max(1, round(Int, n_steps / 15))
    for i in 1:step_stride:n_steps
        state = lambda_fast[i] < -0.1 / p.epsilon ? "Attracting" :
                (lambda_fast[i] > 0.1 / p.epsilon ? "Repelling" : "Fold Trans.")
        @printf("%8.2f | %10.4f | %12.4e | %12.4e | %12.4e | %-15s\n",
            t_arr[i], S_arr[i], e_arr[i], lambda_fast[i], manifold_residual[i], state)
    end
    println("="^80)

    return S_fold_reg, e_fold_reg, S_arr, e_arr, lambda_fast, manifold_residual
end

function main()
    simulate_sbl_dynamics()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end