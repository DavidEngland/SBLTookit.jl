# ==============================================================================
# SBL Vertical Profile Processing Pipeline (Stages 1-4)
# Implements Log-Coordinate Tikhonov-Morozov Spline Smoothing,
# Downstream Businger-Dyer Richardson-Number Inversion, and
# Stage 4 Uncertainty-Aware Iteratively Reweighted Least Squares (IRLS).
# ==============================================================================

using LinearAlgebra
using Printf

# ==============================================================================
# Core Spline Solver and Morozov Bisection
# ==============================================================================

"""
    solve_smoothing_spline(xi, y, sigma, alpha)

Solves a Tikhonov-regularized cubic smoothing spline on the log-height grid `xi`.
Returns the smoothed values, analytical derivatives, and the propagated covariance matrix.
"""
function solve_smoothing_spline(xi::Vector{Float64}, y::Vector{Float64}, sigma::Vector{Float64}, alpha::Float64)
    N = length(xi)
    h = diff(xi)
    
    # 1. Construct the second-difference operator Q (N x N-2)
    Q = zeros(N, N-2)
    for j in 1:(N-2)
        Q[j, j] = 1.0 / h[j]
        Q[j+1, j] = -(1.0 / h[j] + 1.0 / h[j+1])
        Q[j+2, j] = 1.0 / h[j+1]
    end
    
    # 2. Construct the tridiagonal second-derivative continuity matrix R (N-2 x N-2)
    R = zeros(N-2, N-2)
    for j in 1:(N-2)
        R[j, j] = (h[j] + h[j+1]) / 3.0
        if j < N-2
            R[j, j+1] = h[j+1] / 6.0
            R[j+1, j] = h[j+1] / 6.0
        end
    end
    
    # 3. Formulate the Tikhonov linear system with weight diagonal W
    W = Diagonal(1.0 ./ (sigma .^ 2))
    R_inv_QT = R \ Q'
    K = Q * R_inv_QT  # Roughness penalty matrix
    
    A = W + alpha * K
    s = A \ (W * y)  # Smoothed profile values at knots
    
    # 4. Propagate covariance: Sigma_s = A^-1 * W * A^-1
    # Memory-Efficient Covariance Extraction Note:
    # For large systems (N > 1000), explicit dense matrix inversion inv(A)
    # should be replaced with a band-solver factorization or column-by-column
    # solve (A \ e_i) to extract the diagonal covariance elements.
    # Since N is small here (N = 10), dense inv(A) is highly efficient and exact.
    A_inv = inv(A)
    Sigma_s = A_inv * W * A_inv
    
    # 5. Build G_mat such that g = G_mat * s (second derivatives at knots)
    G_mat = zeros(N, N)
    G_mat[2:N-1, :] = R_inv_QT
    
    # 6. Construct derivative matrix D (N x N) for analytical s' = D * s
    D = zeros(N, N)
    for i in 1:(N-1)
        D[i, i] = -1.0 / h[i]
        D[i, i+1] = 1.0 / h[i]
        D[i, :] -= (h[i] / 3.0) * G_mat[i, :] + (h[i] / 6.0) * G_mat[i+1, :]
    end
    D[N, N-1] = -1.0 / h[N-1]
    D[N, N] = 1.0 / h[N-1]
    D[N, :] += (h[N-1] / 6.0) * G_mat[N-1, :] + (h[N-1] / 3.0) * G_mat[N, : ]
    
    s_prime = D * s
    Sigma_s_prime = D * Sigma_s * D'
    
    return s, s_prime, Sigma_s_prime
end

"""
    fit_with_morozov(xi, y, sigma)

Finds the optimal regularization parameter alpha using Morozov's Discrepancy Principle
via a robust, self-expanding bisection search. Targets a chi-squared discrepancy equal to N.
"""
function fit_with_morozov(xi::Vector{Float64}, y::Vector{Float64}, sigma::Vector{Float64})
    N = length(xi)
    low = -15.0
    high = 15.0
    tol = 1e-4
    max_iter = 100
    
    # Helper to calculate chi-squared discrepancy at a given exponent
    calc_chi2(exp_val) = sum(((solve_smoothing_spline(xi, y, sigma, 10.0^exp_val)[1] .- y) ./ sigma) .^ 2)
    
    chi2_low = calc_chi2(low)
    chi2_high = calc_chi2(high)
    
    # Robust self-expanding interval expansion
    while chi2_low > N && low > -30.0
        low -= 5.0
        chi2_low = calc_chi2(low)
    end
    while chi2_high < N && high < 30.0
        high += 5.0
        chi2_high = calc_chi2(high)
    end
    
    if (chi2_low - N) * (chi2_high - N) > 0
        @warn "Morozov target discrepancy N=$N is not bracketed by [10^$low, 10^$high]. Using closest boundary."
        alpha = abs(chi2_low - N) < abs(chi2_high - N) ? 10.0^low : 10.0^high
        s, s_prime, Sigma_s_prime = solve_smoothing_spline(xi, y, sigma, alpha)
        return s, s_prime, Sigma_s_prime, alpha
    end
    
    local s, s_prime, Sigma_s_prime, alpha
    for iter in 1:max_iter
        mid = (low + high) / 2.0
        alpha = 10.0^mid
        s, s_prime, Sigma_s_prime = solve_smoothing_spline(xi, y, sigma, alpha)
        chi2 = sum(((s .- y) ./ sigma) .^ 2)
        
        if abs(chi2 - N) < tol
            break
        elseif chi2 > N
            high = mid  # Over-smoothed: reduce alpha
        else
            low = mid   # Under-smoothed: increase alpha
        end
    end
    return s, s_prime, Sigma_s_prime, alpha
end

# ==============================================================================
# SBL Pipeline Execution and Stage 4 IRLS Implementation
# ==============================================================================

struct PipelineOutput
    z::Vector{Float64}
    
    # Baseline (Stages 1-3)
    theta_smooth_base::Vector{Float64}
    U_smooth_base::Vector{Float64}
    theta_z_base::Vector{Float64}
    U_z_base::Vector{Float64}
    Ri_g_base::Vector{Float64}
    zeta_base::Vector{Float64}
    sigma_zeta_base::Vector{Float64}
    alpha_theta_base::Float64
    alpha_U_base::Float64
    
    # IRLS Refined (Stage 4)
    theta_smooth_irls::Vector{Float64}
    U_smooth_irls::Vector{Float64}
    theta_z_irls::Vector{Float64}
    U_z_irls::Vector{Float64}
    Ri_g_irls::Vector{Float64}
    zeta_irls::Vector{Float64}
    sigma_zeta_irls::Vector{Float64}
    alpha_theta_irls::Float64
    alpha_U_irls::Float64
    w_final::Vector{Float64}
    
    # Diagnostics
    kappa_zeta::Vector{Float64}
    eta_theta::Float64
    eta_U::Float64
    is_instability_flagged::Bool
    irls_iterations::Int
    irls_converged::Bool
end

"""
    run_sbl_pipeline(z, theta_raw, U_raw, delta_theta, delta_U; ...)

Runs the complete 4-Stage SBL Profile processing engine.
* Stage 1-3: Baseline log-coordinate smoothing spline with Businger-Dyer inversion.
* Stage 4: Uncertainty-Aware Damped Iteratively Reweighted Least Squares (IRLS).
"""
function run_sbl_pipeline(\n    z::Vector{Float64}, \n    theta_raw::Vector{Float64}, \n    U_raw::Vector{Float64}, \n    delta_theta::Float64, \n    delta_U::Float64; \n    g=9.81, \n    theta_ref=290.0, \n    Ri_guard=0.19, \n    sigma_zeta_ref=0.05, \n    gamma=0.3, \n    max_irls_iter=30, \n    tol_irls=1e-4\n)
    N = length(z)
    z_0 = 0.1 # reference surface roughness scale
    xi = log.(z ./ z_0)
    
    # ==========================================================================
    # RUN BASELINE (STAGES 1-3)
    # ==========================================================================
    sigma_theta_base = fill(delta_theta, N)
    sigma_U_base = fill(delta_U, N)
    
    theta_smooth_base, theta_xi_prime_base, Sigma_theta_xi_prime_base, alpha_theta_base = fit_with_morozov(xi, theta_raw, sigma_theta_base)
    U_smooth_base, U_xi_prime_base, Sigma_U_xi_prime_base, alpha_U_base = fit_with_morozov(xi, U_raw, sigma_U_base)
    
    theta_z_base = theta_xi_prime_base ./ z
    U_z_base = U_xi_prime_base ./ z
    sigma_theta_z_base = sqrt.(diag(Sigma_theta_xi_prime_base)) ./ z
    sigma_U_z_base = sqrt.(diag(Sigma_U_xi_prime_base)) ./ z
    
    Ri_g_base = (g / theta_ref) .* theta_z_base ./ (U_z_base .^ 2)
    zeta_base = zeros(N)
    sigma_zeta_base = zeros(N)
    
    for i in 1:N
        Ri_val = Ri_g_base[i]
        Ri_guarded = min(Ri_val, Ri_guard)
        
        dzeta_dRi = 0.0
        if Ri_val >= 0.0
            zeta_base[i] = Ri_guarded / (1.0 - 5.0 * Ri_guarded)
            dzeta_dRi = 1.0 / ((1.0 - 5.0 * Ri_guarded)^2)
        else
            zeta_base[i] = Ri_val
            dzeta_dRi = 1.0
        end
        
        U_z_guarded = max(abs(U_z_base[i]), 1e-4) # Low-shear guard
        dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z_guarded^2))
        dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z_base[i] / (theta_ref * U_z_guarded^3))
        
        sigma_zeta_base[i] = sqrt(
            (dzeta_dtheta_z * sigma_theta_z_base[i])^2 + 
            (dzeta_dU_z * sigma_U_z_base[i])^2
        )
    end
    
    # ==========================================================================
    # STAGE 4: UNCERTAINTY-AWARE IRLS LOOP
    # ==========================================================================
    w = ones(N) # Initialize weights
    
    local theta_smooth_irls, U_smooth_irls, theta_z_irls, U_z_irls, Ri_g_irls, zeta_irls, sigma_zeta_irls
    local alpha_theta_irls = alpha_theta_base
    local alpha_U_irls = alpha_U_base
    local irls_iter = 0
    local irls_converged = false
    
    zeta_irls = copy(zeta_base)
    sigma_zeta_irls = copy(sigma_zeta_base)
    
    while irls_iter < max_irls_iter
        irls_iter += 1
        
        # 1. Update data-noise levels using the weights (sigma = delta / sqrt(w))
        sigma_theta_w = delta_theta ./ sqrt.(w)
        sigma_U_w = delta_U ./ sqrt.(w)
        
        # 2. Refit profiles with the new observation weights
        theta_smooth_irls, theta_xi_prime, Sigma_theta_xi_prime, alpha_theta_irls = fit_with_morozov(xi, theta_raw, sigma_theta_w)
        U_smooth_irls, U_xi_prime, Sigma_U_xi_prime, alpha_U_irls = fit_with_morozov(xi, U_raw, sigma_U_w)
        
        # 3. Extract physical derivatives
        theta_z_irls = theta_xi_prime ./ z
        U_z_irls = U_xi_prime ./ z
        sigma_theta_z_irls = sqrt.(diag(Sigma_theta_xi_prime)) ./ z
        sigma_U_z_irls = sqrt.(diag(Sigma_U_xi_prime)) ./ z
        
        # 4. Map Richardson number and non-dimensional stability coordinate
        Ri_g_irls = (g / theta_ref) .* theta_z_irls ./ (U_z_irls .^ 2)
        zeta_new = zeros(N)
        sigma_zeta_new = zeros(N)
        
        for i in 1:N
            Ri_val = Ri_g_irls[i]
            Ri_guarded = min(Ri_val, Ri_guard)
            
            dzeta_dRi = 0.0
            if Ri_val >= 0.0
                zeta_new[i] = Ri_guarded / (1.0 - 5.0 * Ri_guarded)
                dzeta_dRi = 1.0 / ((1.0 - 5.0 * Ri_guarded)^2)
            else
                zeta_new[i] = Ri_val
                dzeta_dRi = 1.0
            end
            
            U_z_guarded = max(abs(U_z_irls[i]), 1e-4) # Low-shear guard
            dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z_guarded^2))
            dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z_irls[i] / (theta_ref * U_z_guarded^3))
            
            sigma_zeta_new[i] = sqrt(
                (dzeta_dtheta_z * sigma_theta_z_irls[i])^2 + 
                (dzeta_dU_z * sigma_U_z_irls[i])^2
            )
        end
        
        # 5. Compute new target weights based on propagated uncertainty
        w_calc = zeros(N)
        for i in 1:N
            w_calc[i] = 1.0 / (1.0 + (sigma_zeta_new[i] / sigma_zeta_ref)^2)
            w_calc[i] = max(w_calc[i], 1e-4) # Floor on data weights
        end
        
        # 6. Apply exponential damping (relaxation) to weight updates
        w_new = (1.0 - gamma) .* w .+ gamma .* w_calc
        
        # 7. Evaluate dual-convergence criteria on both state and weights
        norm_zeta = norm(zeta_irls, Inf)
        norm_w = norm(w, Inf)
        
        diff_zeta = norm(zeta_new .- zeta_irls, Inf) / (1.0 + norm_zeta)
        diff_w = norm(w_new .- w, Inf) / (1.0 + norm_w)
        
        # Update references for the next iteration
        zeta_irls = zeta_new
        sigma_zeta_irls = sigma_zeta_new
        w_prev = copy(w)
        w = w_new
        
        if diff_zeta < tol_irls && diff_w < tol_irls
            irls_converged = true
            break
        end
    end
    
    # Final diagnostic calculations
    Ri_guarded_base = min.(Ri_g_base, Ri_guard)
    kappa_zeta = 1.0 ./ ((1.0 .- 5.0 .* Ri_guarded_base) .^ 2)
    
    eta_theta = alpha_theta_irls / alpha_theta_base
    eta_U = alpha_U_irls / alpha_U_base
    is_instability_flagged = (eta_theta > 3.0 || eta_U > 3.0)
    
    return PipelineOutput(
        z,
        theta_smooth_base, U_smooth_base, theta_z_base, U_z_base, Ri_g_base, zeta_base, sigma_zeta_base, alpha_theta_base, alpha_U_base,
        theta_smooth_irls, U_smooth_irls, theta_z_irls, U_z_irls, Ri_g_irls, zeta_irls, sigma_zeta_irls, alpha_theta_irls, alpha_U_irls, w,
        kappa_zeta, eta_theta, eta_U, is_instability_flagged, irls_iter, irls_converged
    )
end

# ==============================================================================
# Numerical Verification Block
# ==============================================================================

# Define SBL met tower grid
z = [0.5, 1.0, 1.8, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0, 14.0]

# True physical profiles
theta_true = 285.0 .+ 2.0 .* log.(z ./ 0.1) .- 0.05 .* z
U_true = 0.4 .* log.(z ./ 0.1) .+ 0.1 .* z

# Noise parameters
delta_theta = 0.05 # fine-wire thermocouple (K)
delta_U = 0.02     # sonic anemometer (m/s)

# Generate synthetic observations with Gaussian noise
using Random
Random.seed!(42)
theta_obs = theta_true .+ delta_theta .* randn(length(z))
U_obs = U_true .+ delta_U .* randn(length(z))

# Execute SBL processing engine (Stages 1-4)
res = run_sbl_pipeline(z, theta_obs, U_obs, delta_theta, delta_U; sigma_zeta_ref=0.08)

# Format and report results
println("=====================================================================================")
println("SBL STAGES 1-4 COMPREHENSIVE PIPELINE AUDIT")
println("=====================================================================================")
@printf("IRLS Convergence:   %-15s | Iterations: %d\n", res.irls_converged ? "CONVERGED" : "FAILED", res.irls_iterations)
@printf("Alpha Theta (Base): %-15.4e | Alpha Theta (IRLS): %.4e (Ratio: %.2f)\n", res.alpha_theta_base, res.alpha_theta_irls, res.eta_theta)
@printf("Alpha Wind (Base):  %-15.4e | Alpha Wind (IRLS):  %.4e (Ratio: %.2f)\n", res.alpha_U_base, res.alpha_U_irls, res.eta_U)
@printf("Diagnostic Flag:    %-15s\n", res.is_instability_flagged ? "FLAGGED (Instability Suspect)" : "NORMAL")
println("=====================================================================================")
@printf("%-6s | %-12s | %-12s | %-10s | %-10s | %-10s\n", 
        "z (m)", "w_i (Weight)", "Ri_g (Base)", "Ri_g (IRLS)", "zeta_base", "zeta_irls")
println("-------------------------------------------------------------------------------------")
for i in 1:length(z)
    @printf("%-6.1f | %-12.6f | %-12.4f | %-10.4f | %-10.4f | %-10.4f\n",
            res.z[i], res.w_final[i], res.Ri_g_base[i], res.Ri_g_irls[i], res.zeta_base[i], res.zeta_irls[i])
end
println("=====================================================================================")
