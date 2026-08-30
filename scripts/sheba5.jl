# ==============================================================================
# STAGE 5: LATIN HYPERCUBE SAMPLING (LHS) UQ ENGINE
# ==============================================================================

using LatinHypercubeSampling

"""
    run_lhs_uq(z, theta_obs, U_obs, delta_theta, delta_U; N_samples=200)

Draws N_samples from a joint Gaussian noise floor around the raw profiles, passes
each full realization through Stages 1-4, and extracts percentile bounds (P16, P50, P84).
"""
function run_lhs_uq(
    z::Vector{Float64},
    theta_obs::Vector{Float64},
    U_obs::Vector{Float64},
    delta_theta::Float64,
    delta_U::Float64;
    N_samples=200
)
    N_levels = length(z)
    
    # Generate an LHS plan scaled to standard normal distributions N(0,1)
    plan, _ = LHCoptim(N_samples, 2 * N_levels, 100)
    # Map uniform hypercube [0,1] to standard normal quantiles
    using Distributions
    norm_dist = Normal(0.0, 1.0)
    Z_scores = quantile.(norm_dist, plan)
    
    zeta_ensemble = zeros(N_levels, N_samples)
    
    for s in 1:N_samples
        # Perturb raw observations with independent Gaussian noise
        noise_theta = Z_scores[s, 1:N_levels] .* delta_theta
        noise_U = Z_scores[s, (N_levels+1):(2*N_levels)] .* delta_U

        theta_pert = theta_obs .+ noise_theta
        U_pert = U_obs .+ noise_U

        # Run full 4-stage pipeline on perturbed profile
        res_s = run_sbl_pipeline(z, theta_pert, U_pert, delta_theta, delta_U)
        zeta_ensemble[:, s] = res_s.zeta_irls
    end
    
    # Compute non-parametric percentiles across the ensemble
    P16 = [quantile(zeta_ensemble[i, :], 0.16) for i in 1:N_levels]
    P50 = [quantile(zeta_ensemble[i, :], 0.50) for i in 1:N_levels]
    P84 = [quantile(zeta_ensemble[i, :], 0.84) for i in 1:N_levels]
    
    return P16, P50, P84
end

# ==============================================================================
# STAGE 6: GRACHEV ET AL. (2007) SHEBA CLOSURE MODULE
# ==============================================================================

struct GrachevClosure end

"""
    eval_grachev_phi(zeta)

Computes the Grachev et al. (2007) stability functions phi_m and phi_h for stable conditions.
"""
function eval_grachev_phi(zeta::Float64)
    a_m, b_m = 5.0, 5.0 / 6.5
    a_h, b_h, c_h = 5.0, 5.0, 3.0
    
    phi_m = 1.0 + a_m * (zeta * (1.0 + zeta)^(1.0/3.0)) / (1.0 + b_m * zeta)
    phi_h = 1.0 + (a_h * zeta + b_h * zeta^2) / (1.0 + c_h * zeta + zeta^2)
    
    return phi_m, phi_h
end

"""
    invert_grachev_Ri(Ri_g; tol=1e-6, max_iter=50)

Inverts Ri_g to zeta for the Grachev (2007) closure via Newton-Raphson root finding on:
F(zeta) = zeta * phi_h(zeta) / (phi_m(zeta)^2) - Ri_g = 0
Returns both zeta and the analytical derivative dzeta/dRi_g = 1 / R'(zeta).
"""
function invert_grachev_Ri(Ri_g::Float64; tol=1e-6, max_iter=50)
    if Ri_g < 0.0
        return Ri_g, 1.0 # Unstable branch fallback
    end
    
    zeta = Ri_g # Initial guess
    
    for iter in 1:max_iter
        phi_m, phi_h = eval_grachev_phi(zeta)
        R_val = zeta * phi_h / (phi_m^2)
        
        f = R_val - Ri_g
        if abs(f) < tol
            break
        end
        
        # Numerical/AD derivative of R(zeta) for exact step and dzeta/dRi_g
        dzeta_step = 1e-5
        phi_m_p, phi_h_p = eval_grachev_phi(zeta + dzeta_step)
        R_val_p = (zeta + dzeta_step) * phi_h_p / (phi_m_p^2)
        
        dR_dzeta = (R_val_p - R_val) / dzeta_step
        zeta -= f / dR_dzeta
        zeta = max(zeta, 0.0)
    end
    
    # Final derivative evaluation for sensitivity matrix
    phi_m, phi_h = eval_grachev_phi(zeta)
    dzeta_step = 1e-5
    phi_m_p, phi_h_p = eval_grachev_phi(zeta + dzeta_step)
    R_val_p = (zeta + dzeta_step) * phi_h_p / (phi_m_p^2)
    R_val = zeta * phi_h / (phi_m^2)
    
    dR_dzeta = (R_val_p - R_val) / dzeta_step
    dzeta_dRi = 1.0 / max(dR_dzeta, 1e-6)
    
    return zeta, dzeta_dRi
end