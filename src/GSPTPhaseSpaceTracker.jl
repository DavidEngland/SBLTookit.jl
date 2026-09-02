module GSPTPhaseSpaceTracker

using Statistics
using LinearAlgebra
using Printf

export CASES99SiteParams, GSPTState, evaluate_fast_manifold!, step_stiffness_aware!, run_phase_space_tracking

"""
    CASES99SiteParams

Type-stable configuration parameter struct for the CASES-99 Night 991018 site parameters.
Can be overridden with user-specified values at runtime.
"""
Base.@kwdef struct CASES99SiteParams
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
    dt_base::Float64 = 0.001       # Base timestep size (s)
    safety_eta::Float64 = 0.50     # Timestep safety factor for stiffness sub-stepping
end

"""
    GSPTState

Mutable struct representing the in-place state of the fast-slow SBL trajectory tracker.
Ensures zero heap allocations inside the hot integration loop.
"""
mutable struct GSPTState
    e::Float64                     # Turbulent Kinetic Energy (fast variable)
    S::Float64                     # Wind shear (slow variable)
    t::Float64                     # Cumulative time
    F::Float64                     # Fast subsystem vector field value f(e, S)
    g::Float64                     # Slow subsystem vector field value g(e, S)
    Fe::Float64                    # Fast subsystem Jacobian derivative w.r.t e (lambda_fast)
    FS::Float64                    # Fast subsystem Jacobian derivative w.r.t S
    Fee::Float64                   # 2nd derivative w.r.t e
    FeS::Float64                   # Cross derivative w.r.t e and S
    d_perp::Float64                # Transverse geometric distance to critical manifold
    dist_fold::Float64             # Geometric distance to the saddle-node fold boundary
end

"""
    evaluate_fast_manifold!(state::GSPTState, p::CASES99SiteParams)

Computes the fast vector field, slow vector field, exact local fast eigenvalue (Fe),
and analytical gradient components (FS, Fee, FeS) in-place with zero heap allocations.
"""
function evaluate_fast_manifold!(state::GSPTState, p::CASES99SiteParams)
    e = state.e
    S = state.S
    Seff = sqrt(S^2 + p.eta^2)

    # Buoyancy flux and derivatives
    denom_e = e^2 + p.delta_reg^2
    B = p.B0_max * (e^2) / denom_e
    Be = 2.0 * p.B0_max * e * (p.delta_reg^2) / (denom_e^2)
    Bee = 2.0 * p.B0_max * (p.delta_reg^2) * (p.delta_reg^2 - 3.0 * e^2) / (denom_e^3)

    # Dissipation and derivatives
    denom_S = Seff^2 + p.a
    Diss = (e^3) * (Seff^2) / (p.l0 * denom_S)

    # Update state vector fields
    state.F = p.l0 * e * (Seff^2) - B - Diss
    state.g = p.G0 - p.gamma_s * e * S - p.r_s * S

    # Analytical Jacobian derivatives
    state.Fe = p.l0 * (Seff^2) - Be - 3.0 * (e^2) * (Seff^2) / (p.l0 * denom_S)
    state.Fee = - Bee - 6.0 * e * (Seff^2) / (p.l0 * denom_S)

    # Derivatives w.r.t S (exact chain-rule treatment)
    dS_eff_dS = S / Seff
    dF_dSeff = 2.0 * p.l0 * e * Seff - (2.0 * Seff * p.a * (e^3)) / (p.l0 * (denom_S^2))
    dFe_dSeff = 2.0 * p.l0 * Seff - (6.0 * Seff * p.a * (e^2)) / (p.l0 * (denom_S^2))

    state.FS = dF_dSeff * dS_eff_dS
    state.FeS = dFe_dSeff * dS_eff_dS

    # Calculate geometric transverse manifold distance
    state.d_perp = abs(state.F) / (sqrt(state.Fe^2 + state.FS^2) + 1e-12)
end

"""
    step_stiffness_aware!(state::GSPTState, p::CASES99SiteParams, S_fold::Float64, e_fold::Float64)

Performs a single in-place RK4 integration step utilizing a stiffness-aware sub-stepping bound
to prevent unphysical overshoot into negative TKE space during fast relaxation phases.
"""
function step_stiffness_aware!(state::GSPTState, p::CASES99SiteParams, S_fold::Float64, e_fold::Float64)
    # Evaluate current manifold metrics
    evaluate_fast_manifold!(state, p)

    # Stiffness-aware adaptive sub-stepping
    lambda_fast = state.Fe / p.epsilon
    dt_adaptive = min(p.dt_base, p.safety_eta / (abs(lambda_fast) + 1e-12))

    # In-place RK4 integration step
    # k1
    e1, S1 = state.e, state.S
    F1 = state.F
    g1 = state.g
    de1 = F1 / p.epsilon
    dS1 = g1

    # k2
    e2 = max(e1 + 0.5 * dt_adaptive * de1, 1e-6)
    S2 = max(S1 + 0.5 * dt_adaptive * dS1, 0.0)

    state.e = e2
    state.S = S2
    evaluate_fast_manifold!(state, p)
    F2 = state.F
    g2 = state.g
    de2 = F2 / p.epsilon
    dS2 = g2

    # k3
    e3 = max(e1 + 0.5 * dt_adaptive * de2, 1e-6)
    S3 = max(S1 + 0.5 * dt_adaptive * dS2, 0.0)

    state.e = e3
    state.S = S3
    evaluate_fast_manifold!(state, p)
    F3 = state.F
    g3 = state.g
    de3 = F3 / p.epsilon
    dS3 = g3

    # k4
    e4 = max(e1 + dt_adaptive * de3, 1e-6)
    S4 = max(S1 + dt_adaptive * dS3, 0.0)

    state.e = e4
    state.S = S4
    evaluate_fast_manifold!(state, p)
    F4 = state.F
    g4 = state.g
    de4 = F4 / p.epsilon
    dS4 = g4

    # Final in-place updates
    state.e = max(e1 + (dt_adaptive / 6.0) * (de1 + 2.0*de2 + 2.0*de3 + de4), 1e-6)
    state.S = max(S1 + (dt_adaptive / 6.0) * (dS1 + 2.0*dS2 + 2.0*dS3 + dS4), 0.0)
    state.t += dt_adaptive

    # Re-evaluate final manifold metrics
    evaluate_fast_manifold!(state, p)

    # Update distance to fold boundary
    state.dist_fold = sqrt((state.S - S_fold)^2 + (state.e - e_fold)^2)
end

"""
    run_phase_space_tracking(p::CASES99SiteParams, S_fold::Float64, e_fold::Float64; t_max=80.0)

Executes the complete trajectory tracking run and exports the high-precision state-space coordinates.
"""
function run_phase_space_tracking(p::CASES99SiteParams, S_fold::Float64, e_fold::Float64; t_max=80.0)
    # Initialize on manifold (Fe < 0)
    e0 = e_fold * 1.5  # Seeding guess
    state = GSPTState(e0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Newton search for on-manifold initialization
    for i in 1:100
        evaluate_fast_manifold!(state, p)
        if abs(state.F) < 1e-12
            break
        end
        state.e = max(state.e - state.F / state.Fe, 1e-6)
    end

    # Pre-allocate timeseries arrays
    max_steps = round(Int, t_max / p.dt_base) * 2
    t_hist = zeros(max_steps)
    e_hist = zeros(max_steps)
    S_hist = zeros(max_steps)
    Fe_hist = zeros(max_steps)
    FS_hist = zeros(max_steps)
    d_perp_hist = zeros(max_steps)
    dist_fold_hist = zeros(max_steps)

    step_idx = 1
    t_hist[1] = state.t
    e_hist[1] = state.e
    S_hist[1] = state.S
    Fe_hist[1] = state.Fe
    FS_hist[1] = state.FS
    d_perp_hist[1] = state.d_perp
    dist_fold_hist[1] = state.dist_fold

    # Simulation Loop
    while state.t < t_max
        step_idx += 1
        step_stiffness_aware!(state, p, S_fold, e_fold)

        t_hist[step_idx] = state.t
        e_hist[step_idx] = state.e
        S_hist[step_idx] = state.S
        Fe_hist[step_idx] = state.Fe
        FS_hist[step_idx] = state.FS
        d_perp_hist[step_idx] = state.d_perp
        dist_fold_hist[step_idx] = state.dist_fold
    end

    # Slice to actual integrated size
    t_out = t_hist[1:step_idx]
    e_out = e_hist[1:step_idx]
    S_out = S_hist[1:step_idx]
    Fe_out = Fe_hist[1:step_idx]
    FS_out = FS_hist[1:step_idx]
    d_perp_out = d_perp_hist[1:step_idx]
    dist_fold_out = dist_fold_hist[1:step_idx]

    return t_out, e_out, S_out, Fe_out, FS_out, d_perp_out, dist_fold_out
end

end # module
