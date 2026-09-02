# SBLToolKit.jl: GSPT Saddle-Node Trajectory Simulation
# src/GSPTSaddleNodeTrajectory.jl
module GSPTSaddleNodeTrajectory

using Statistics
using LinearAlgebra
using Printf

export SBLParams, fast_manifold_all, solve_asymptotic_fold, solve_regularized_fold, initialize_on_manifold, simulate_sbl_dynamics

"""
    SBLParams

Type-stable configuration parameter struct for SBL fast-slow saddle-node trajectory simulation.
"""
Base.@kwdef struct SBLParams
    l0::Float64 = 1.0              # Mixing length scale (m)
    beta::Float64 = 5.0            # Monin-Obukhov stable similarity constant
    N2::Float64 = 0.10             # Buoyancy frequency (s^-2)
    a::Float64 = 0.5               # beta * N2 parameter
    B0_max::Float64 = 0.05         # Maximum buoyant destruction limit (m^2/s^3)
    delta_reg::Float64 = 0.01      # Buoyancy regularization floor (m^2/s^2)
    eta::Float64 = 1e-4            # Shear coordinate smoothing parameter (s^-1)
    G0::Float64 = 0.3              # Geostrophic wind shear forcing (s^-2)
    gamma_s::Float64 = 1.8         # Turbulent damping feedback (s*m^-2)
    r_s::Float64 = 0.15            # Linear shear decay rate (s^-1)
    epsilon::Float64 = 0.05        # Timescale separation parameter
end

# Copy constructor for easy delta-sweeps
function SBLParams(p::SBLParams; delta_reg=p.delta_reg)
    return SBLParams(
        l0=p.l0, beta=p.beta, N2=p.N2, a=p.a, B0_max=p.B0_max,
        delta_reg=delta_reg, eta=p.eta, G0=p.G0, gamma_s=p.gamma_s,
        r_s=p.r_s, epsilon=p.epsilon
    )
end

"""
    fast_manifold_all(e, S, p::SBLParams)

Computes fast subsystem vector field F(e, S) and all first and second order analytical
derivatives: F_e, F_S, F_ee, F_eS. Integrates smooth S_eff coordinate regularization.
"""
function fast_manifold_all(e::Float64, S::Float64, p::SBLParams)
    Seff = sqrt(S^2 + p.eta^2)

    # Buoyancy flux and derivatives
    denom_e = e^2 + p.delta_reg^2
    B = p.B0_max * (e^2) / denom_e
    Be = 2.0 * p.B0_max * e * (p.delta_reg^2) / (denom_e^2)
    Bee = 2.0 * p.B0_max * (p.delta_reg^2) * (p.delta_reg^2 - 3.0 * e^2) / (denom_e^3)

    # Dissipation and derivatives
    denom_S = Seff^2 + p.a
    Diss = (e^3) * (Seff^2) / (p.l0 * denom_S)

    F = p.l0 * e * (Seff^2) - B - Diss
    Fe = p.l0 * (Seff^2) - Be - 3.0 * (e^2) * (Seff^2) / (p.l0 * denom_S)
    Fee = - Bee - 6.0 * e * (Seff^2) / (p.l0 * denom_S)

    # Derivatives of S_eff with respect to S (Chain rule)
    FS = (2.0 * p.l0 * e - 2.0 * p.a * (e^3) / (p.l0 * (denom_S^2))) * S
    FeS = (2.0 * p.l0 - 6.0 * p.a * (e^2) / (p.l0 * (denom_S^2))) * S

    return F, Fe, FS, Fee, FeS
end

"""
    solve_asymptotic_fold(p::SBLParams)

Analytical closed-form fold solver for the unregularized (delta_reg = 0, eta = 0) system
using Cardano's formula.
"""
function solve_asymptotic_fold(p::SBLParams)
    a = p.a
    C = 27.0 * (p.B0_max^2) / (4.0 * (p.l0^4))

    p_c = - (a^2) / 3.0
    q_c = (2.0 * (a^3)) / 27.0 - C

    D = (q_c^2) / 4.0 + (p_c^3) / 27.0
    if D < 0
        phi = acos(clamp(-q_c / (2.0 * sqrt(-(p_c^3) / 27.0)), -1.0, 1.0))
        r_c = 2.0 * sqrt(-p_c / 3.0)
        y = r_c * cos(phi / 3.0)
        x_fold = y - a / 3.0
    else
        u_c = cbrt(-q_c / 2.0 + sqrt(D))
        v_c = cbrt(-q_c / 2.0 - sqrt(D))
        x_fold = u_c + v_c - a / 3.0
    end

    S_fold = sqrt(x_fold)
    e_fold = (p.l0 / sqrt(3.0)) * sqrt(S_fold^2 + a)
    return S_fold, e_fold
end

"""
    solve_regularized_fold(p::SBLParams; e_guess=0.47, S_guess=0.40, max_iter=100, tol=1e-12)

2D Newton-Raphson solver to find the precise regularized fold point satisfying F = 0, Fe = 0.
"""
function solve_regularized_fold(p::SBLParams; e_guess=0.47, S_guess=0.40, max_iter=100, tol=1e-12)
    e = e_guess
    S = S_guess
    for i in 1:max_iter
        F, Fe, FS, Fee, FeS = fast_manifold_all(e, S, p)
        R = [F, Fe]
        if norm(R) < tol
            return S, e, true
        end
        J = [Fe FS; Fee FeS]
        try
            delta = J \ -R
            e += delta[1]
            S += delta[2]
        catch
            return NaN, np.nan, false
        end
    end
    return S, e, false
end

"""
    initialize_on_manifold(S_0, p::SBLParams; e_guess=0.8, max_iter=100, tol=1e-12)

Newton solver to place initial state directly on the attracting slow manifold (Fe < 0).
"""
function initialize_on_manifold(S_0::Float64, p::SBLParams; e_guess=0.8, max_iter=100, tol=1e-12)
    e = e_guess
    for i in 1:max_iter
        F, Fe, _, _, _ = fast_manifold_all(e, S_0, p)
        if abs(F) < tol
            if Fe < 0
                return e, true
            else
                e = e_guess + 0.1
                continue
            end
        end
        e -= F / Fe
    end
    return e, false
end

"""
    simulate_sbl_dynamics(p::SBLParams; t_max=80.0, dt=0.001)

Performs high-fidelity RK4 slow-fast integration and outputs flat trajectory arrays alongside
the exact indices of fold crossings, hyperbolicity loss, and dynamic collapses.
"""
function simulate_sbl_dynamics(p::SBLParams; t_max=80.0, dt=0.001)
    # Solve fold locations
    S_fold, e_fold, _ = solve_regularized_fold(p)

    # Initialize on manifold
    e0, ok = initialize_on_manifold(1.0, p)
    if !ok
        error("Failed to initialize state on slow manifold.")
    end

    n_steps = round(Int, t_max / dt)
    t_arr = collect(range(0.0, t_max, length=n_steps))
    e_arr = zeros(n_steps)
    S_arr = zeros(n_steps)

    e_arr[1] = e0
    S_arr[1] = 1.0

    # Event lists
    fold_crossings = Tuple{Float64,Float64,Float64}[]
    hyper_losses = Tuple{Float64,Float64,Float64}[]
    dynamic_collapses = Tuple{Float64,Float64,Float64}[]

    for i in 2:n_steps
        t_curr = t_arr[i-1]

        # RK4 substeps
        # k1
        e1, S_1 = e_arr[i-1], S_arr[i-1]
        F1, Fe1, _, _, _ = fast_manifold_all(e1, S_1, p)
        g1 = p.G0 - p.gamma_s * e1 * S_1 - p.r_s * S_1
        de1 = F1 / p.epsilon
        dS1 = g1

        # k2
        e2 = max(e1 + 0.5 * dt * de1, 1e-4)
        S2 = max(S_1 + 0.5 * dt * dS1, 0.0)
        F2, _, _, _, _ = fast_manifold_all(e2, S2, p)
        g2 = p.G0 - p.gamma_s * e2 * S2 - p.r_s * S2
        de2 = F2 / p.epsilon
        dS2 = g2

        # k3
        e3 = max(e1 + 0.5 * dt * de2, 1e-4)
        S3 = max(S_1 + 0.5 * dt * dS2, 0.0)
        F3, _, _, _, _ = fast_manifold_all(e3, S3, p)
        g3 = p.G0 - p.gamma_s * e3 * S3 - p.r_s * S3
        de3 = F3 / p.epsilon
        dS3 = g3

        # k4
        e4 = max(e1 + dt * de3, 1e-4)
        S4 = max(S_1 + dt * dS3, 0.0)
        F4, _, _, _, _ = fast_manifold_all(e4, S4, p)
        g4 = p.G0 - p.gamma_s * e4 * S4 - p.r_s * S4
        de4 = F4 / p.epsilon
        dS4 = g4

        # Update
        e_next = e1 + (dt / 6.0) * (de1 + 2.0*de2 + 2.0*de3 + de4)
        S_next = S_1 + (dt / 6.0) * (dS1 + 2.0*dS2 + 2.0*dS3 + dS4)

        e_arr[i] = max(e_next, 1e-4)
        S_arr[i] = max(S_next, 0.0)

        # Events detection
        if (S_arr[i-1] >= S_fold && S_arr[i] < S_fold) || (S_arr[i-1] <= S_fold && S_arr[i] > S_fold)
            push!(fold_crossings, (t_curr, S_arr[i], e_arr[i]))
        end

        _, Fe_prev, _, _, _ = fast_manifold_all(e_arr[i-1], S_arr[i-1], p)
        _, Fe_curr, _, _, _ = fast_manifold_all(e_arr[i], S_arr[i], p)
        if (Fe_prev < 0 && Fe_curr >= 0) || (Fe_prev > 0 && Fe_curr <= 0)
            push!(hyper_losses, (t_curr, S_arr[i], e_arr[i]))
        end

        if e_arr[i-1] > 0.1 && e_arr[i] <= 0.1
            push!(dynamic_collapses, (t_curr, S_arr[i], e_arr[i]))
        end
    end

    return t_arr, e_arr, S_arr, fold_crossings, hyper_losses, dynamic_collapses
end

end # module
