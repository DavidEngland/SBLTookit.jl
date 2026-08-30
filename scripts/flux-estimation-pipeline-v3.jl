# ==============================================================================
# SBL Vertical Profile Processing Pipeline (Stages 1-6)
# Implements Log-Coordinate Tikhonov-Morozov Spline Smoothing,
# Downstream Businger-Dyer and Grachev (2007) Richardson-Number Inversions,
# Stage 4 Uncertainty-Aware Iteratively Reweighted Least Squares (IRLS),
# Stage 5 Newton-Raphson Numerical Solver for SHEBA-calibrated Physics, and
# Stage 6 GSPT Closure Comparison and Curvature Decomposition.
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
    R_inv_QT = R \\ Q'
    K = Q * R_inv_QT  # Roughness penalty matrix
    
    A = W + alpha * K
    s = A \\ (W * y)  # Smoothed profile values at knots
    
    # 4. Propagate covariance: Sigma_s = A^-1 * W * A^-1
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
    D[N, :] += (h[N-1] / 6.0) * G_mat[N-1, :] + (h[N-1] / 3.0) * G_mat[N, :]
    
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
# Stage 5: Newton-Raphson Solver for Grachev et al. (2007) SHEBA Physics
# ==============================================================================

"""
    phi_m_and_prime(zeta)

Computes the SHEBA-calibrated momentum stability function and its analytical first-order
derivative with respect to zeta based on Grachev et al. (2007).
"""
function phi_m_and_prime(zeta::Float64)
    am = 5.0
    bm = 5.0 / 6.5
    if zeta <= 0.0
        return 1.0, 0.0 # Unstable regime default
    end
    pow1_3 = (1.0 + zeta)^(1.0/3.0)
    pow2_3 = (1.0 + zeta)^(2.0/3.0)
    
    val = 1.0 + am * zeta * pow1_3 / (1.0 + bm * zeta)
    
    num = 1.0 + (4.0/3.0)*zeta + (1.0/3.0)*bm*zeta^2
    den = pow2_3 * (1.0 + bm * zeta)^2
    deriv = am * num / den
    return val, deriv
end

"""
    phi_h_and_prime(zeta)

Computes the SHEBA-calibrated heat stability function and its analytical first-order
derivative with respect to zeta based on Grachev et al. (2007).
"""
function phi_h_and_prime(zeta::Float64)
    ah = 5.0
    bh = 5.0
    ch = 3.0
    if zeta <= 0.0
        return 1.0, 0.0 # Unstable regime default
    end
    val = 1.0 + (ah * zeta + bh * zeta^2) / (1.0 + ch * zeta + zeta^2)
    
    num = ah + 2.0*bh*zeta + (bh*ch - ah)*zeta^2
    den = (1.0 + ch * zeta + zeta^2)^2
    deriv = num / den
    return val, deriv
end

"""
    R_and_prime(zeta)

Computes the physical-gradient-to-similarity mapping R(zeta) = zeta * phi_h(zeta) / phi_m^2(zeta)
and its exact analytical first derivative R'(zeta).
"""
function R_and_prime(zeta::Float64)
    if zeta <= 0.0
        return zeta, 1.0
    end
    pm, pm_prime = phi_m_and_prime(zeta)
    ph, ph_prime = phi_h_and_prime(zeta)
    
    val = zeta * ph / (pm^2)
    
    term1 = (ph + zeta * ph_prime) / (pm^2)
    term2 = 2.0 * zeta * ph * pm_prime / (pm^3)
    prime = term1 - term2
    return val, prime
end

"""
    R_second_derivative_approx(zeta)

Approximates the second derivative R''(zeta) via highly accurate centered finite differences
on the analytical first derivative R'(zeta) with a step h = 1e-5.
"""
function R_second_derivative_approx(zeta::Float64)
    if zeta <= 0.0
        return 0.0
    end
    h = 1e-5
    _, prime_plus = R_and_prime(zeta + h)
    _, prime_minus = R_and_prime(max(zeta - h, 0.0))
    return (prime_plus - prime_minus) / (2.0 * h)
end

"""
    invert_grachev(Ri_g; max_iter=100, tol=1e-8)

Performs a robust, analytically driven Newton-Raphson inversion of the non-monotonic
Grachev et al. (2007) stable similarity relation to resolve zeta and R'(zeta).
"""
function invert_grachev(Ri_g::Float64; max_iter=100, tol=1e-8)
    if Ri_g < 0.0
        return Ri_g, 1.0
    end
    
    # Grachev et al. (2007) has no hard critical Richardson cutoff!
    # At high zeta, R(zeta) grows as zeta^(1/3). We initialize Newton's solver dynamically:
    zeta = Ri_g < 0.1 ? Ri_g : (Ri_g < 0.5 ? 2.0 * Ri_g : 10.0)
    
    for iter in 1:max_iter
        R_val, R_prime = R_and_prime(zeta)
        diff_val = R_val - Ri_g
        if abs(diff_val) < tol
            return zeta, R_prime
        end
        # Newton update
        zeta -= diff_val / R_prime
        if zeta < 0.0
            zeta = 1e-6 # maintain positive constraint
        end
    end
    return zeta, R_and_prime(zeta)[2]
end

# ==============================================================================
# Helper Differentiation Splines for Stage 6 GSPT Curvature Calculations
# ==============================================================================

"""
    extract_zeta_derivatives(xi, z, zeta_profile)

Fits a highly accurate natural cubic spline to the vertical profile of diagnosed zeta
to extract its exact, analytical vertical spatial derivatives zeta_z and zeta_zz.
"""
function extract_zeta_derivatives(xi::Vector{Float64}, z::Vector{Float64}, zeta_profile::Vector{Float64})
    N = length(xi)
    sigma_unit = fill(1e-4, N)
    s, s_prime, _, _ = solve_smoothing_spline(xi, zeta_profile, sigma_unit, 1e-8)
    
    h = diff(xi)
    Q = zeros(N, N-2)
    for j in 1:(N-2)
        Q[j, j] = 1.0 / h[j]
        Q[j+1, j] = -(1.0 / h[j] + 1.0 / h[j+1])
        Q[j+2, j] = 1.0 / h[j+1]
    end
    R = zeros(N-2, N-2)
    for j in 1:(N-2)
        R[j, j] = (h[j] + h[j+1]) / 3.0
        if j < N-2
            R[j, j+1] = h[j+1] / 6.0
            R[j+1, j] = h[j+1] / 6.0
        end
    end
    g_interior = R \\ (Q' * s)
    g = zeros(N)
    g[2:N-1] = g_interior
    
    zeta_z = s_prime ./ z
    zeta_zz = (g .- s_prime) ./ (z .^ 2)
    
    return zeta_z, zeta_zz
end

"""
    extract_Ri_second_derivative(xi, z, Ri_profile)

Fits a highly accurate natural cubic spline to the vertical profile of gradient Richardson numbers
to extract its exact vertical second derivative Ri_zz.
"""
function extract_Ri_second_derivative(xi::Vector{Float64}, z::Vector{Float64}, Ri_profile::Vector{Float64})
    N = length(xi)
    sigma_unit = fill(1e-4, N)
    s, s_prime, _, _ = solve_smoothing_spline(xi, Ri_profile, sigma_unit, 1e-8)
    
    h = diff(xi)
    Q = zeros(N, N-2)
    for j in 1:(N-2)
        Q[j, j] = 1.0 / h[j]
        Q[j+1, j] = -(1.0 / h[j] + 1.0 / h[j+1])
        Q[j+2, j] = 1.0 / h[j+1]
    end
    R = zeros(N-2, N-2)
    for j in 1:(N-2)
        R[j, j] = (h[j] + h[j+1]) / 3.0
        if j < N-2
            R[j, j+1] = h[j+1] / 6.0
            R[j+1, j] = h[j+1] / 6.0
        end
    end
    g_interior = R \\ (Q' * s)
    g = zeros(N)
    g[2:N-1] = g_interior
    
    Ri_zz = (g .- s_prime) ./ (z .^ 2)
    return Ri_zz
end

# ==============================================================================
# SBL Pipeline Execution and Stage 6 Closure Comparison
# ==============================================================================

struct ClosureComparisonResult
    z::Vector{Float64}
    
    # State gradients and observed Richardson number
    theta_smooth::Vector{Float64}
    U_smooth::Vector{Float64}
    theta_z::Vector{Float64}
    U_z::Vector{Float64}
    Ri_g::Vector{Float64}
    Ri_zz_obs::Vector{Float64} # Observed spatial second derivative
    
    # Businger-Dyer (BD) Inversion and Curvature Decomposition
    zeta_BD::Vector{Float64}
    sigma_zeta_BD::Vector{Float64}
    kappa_zeta_BD::Vector{Float64}
    zeta_z_BD::Vector{Float64}
    zeta_zz_BD::Vector{Float64}
    C_const_BD::Vector{Float64}
    C_coord_BD::Vector{Float64}
    Ri_zz_GSPT_BD::Vector{Float64}
    delta_Ri_zz_BD::Vector{Float64} # Residual delta_Ri_zz = Ri_zz_obs - (C_const + C_coord)
    
    # Grachev (G07) Inversion and Curvature Decomposition
    zeta_G07::Vector{Float64}
    sigma_zeta_G07::Vector{Float64}
    kappa_zeta_G07::Vector{Float64}
    zeta_z_G07::Vector{Float64}
    zeta_zz_G07::Vector{Float64}
    C_const_G07::Vector{Float64}
    C_coord_G07::Vector{Float64}
    Ri_zz_GSPT_G07::Vector{Float64}
    delta_Ri_zz_G07::Vector{Float64} # Residual delta_Ri_zz = Ri_zz_obs - (C_const + C_coord)
    
    # Diagnostic Metrics and Flags
    is_high_smoothing::Bool
    eta_theta::Float64
    eta_U::Float64
end

"""
    compare_sbl_closures(z, theta_raw, U_raw, delta_theta, delta_U; ...)

Runs the complete 6-Stage SBL Profile processing engine.
Inverts physical profiles into similarity coordinates under both Businger-Dyer and Grachev (2007) closures,
performs exact analytical error propagation, and computes the vertical spatial curvature decomposition
to isolate true physical turbulence closure discrepancies from coordinate projection artifacts.
"""
function compare_sbl_closures(
    z::Vector{Float64}, 
    theta_raw::Vector{Float64}, 
    U_raw::Vector{Float64}, 
    delta_theta::Float64, 
    delta_U::Float64; 
    g=9.81, 
    theta_ref=290.0, 
    Ri_guard=0.19,
    sigma_zeta_ref=0.05,
    gamma=0.3,
    max_irls_iter=30,
    tol_irls=1e-4
)
    N = length(z)
    z_0 = 4.5e-4 # Reference roughness length for SHEBA campaign (Persson et al., 2001)
    xi = log.(z ./ z_0)
    
    # ==========================================================================
    # STAGES 1-4: UNCERTAINTY-AWARE PROFILE REGULARIZATION (IRLS)
    # ==========================================================================
    # Execute standard Stage 4 IRLS log-coordinate spline to obtain smoothed profiles and gradients
    # We first fit the baseline
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
        dzeta_dRi = Ri_val >= 0.0 ? 1.0 / ((1.0 - 5.0 * Ri_guarded)^2) : 1.0
        U_z_guarded = max(abs(U_z_base[i]), 1e-4)
        dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z_guarded^2))
        dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z_base[i] / (theta_ref * U_z_guarded^3))
        sigma_zeta_base[i] = sqrt((dzeta_dtheta_z * sigma_theta_z_base[i])^2 + (dzeta_dU_z * sigma_U_z_base[i])^2)
    end
    
    # Run the IRLS iteration loop
    w = ones(N)
    local theta_smooth, U_smooth, theta_z, U_z, sigma_theta_z, sigma_U_z, Ri_g
    local alpha_theta_irls = alpha_theta_base
    local alpha_U_irls = alpha_U_base
    local irls_iter = 0
    zeta_irls = copy(zeta_base)
    
    while irls_iter < max_irls_iter
        irls_iter += 1
        sigma_theta_w = delta_theta ./ sqrt.(w)
        sigma_U_w = delta_U ./ sqrt.(w)
        
        theta_smooth, theta_xi_prime, Sigma_theta_xi_prime, alpha_theta_irls = fit_with_morozov(xi, theta_raw, sigma_theta_w)
        U_smooth, U_xi_prime, Sigma_U_xi_prime, alpha_U_irls = fit_with_morozov(xi, U_raw, sigma_U_w)
        
        theta_z = theta_xi_prime ./ z
        U_z = U_xi_prime ./ z
        sigma_theta_z = sqrt.(diag(Sigma_theta_xi_prime)) ./ z
        sigma_U_z = sqrt.(diag(Sigma_U_xi_prime)) ./ z
        
        Ri_g = (g / theta_ref) .* theta_z ./ (U_z .^ 2)
        zeta_new = zeros(N)
        sigma_zeta_new = zeros(N)
        
        for i in 1:N
            Ri_val = Ri_g[i]
            Ri_guarded = min(Ri_val, Ri_guard)
            dzeta_dRi = Ri_val >= 0.0 ? 1.0 / ((1.0 - 5.0 * Ri_guarded)^2) : 1.0
            U_z_guarded = max(abs(U_z[i]), 1e-4)
            dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z_guarded^2))
            dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z[i] / (theta_ref * U_z_guarded^3))
            sigma_zeta_new[i] = sqrt((dzeta_dtheta_z * sigma_theta_z[i])^2 + (dzeta_dU_z * sigma_U_z[i])^2)
            zeta_new[i] = Ri_val >= 0.0 ? Ri_guarded / (1.0 - 5.0 * Ri_guarded) : Ri_val
        end
        
        w_calc = [1.0 / (1.0 + (sigma_zeta_new[i] / sigma_zeta_ref)^2) for i in 1:N]
        w_new = (1.0 - gamma) .* w .+ gamma .* max.(w_calc, 1e-4)
        
        diff_zeta = norm(zeta_new .- zeta_irls, Inf) / (1.0 + norm(zeta_irls, Inf))
        diff_w = norm(w_new .- w, Inf) / (1.0 + norm(w, Inf))
        
        zeta_irls = zeta_new
        w = w_new
        
        if diff_zeta < tol_irls && diff_w < tol_irls
            break
        end
    end
    
    # Compute physical-space second derivative of Gradient Richardson Number: Ri_zz
    Ri_zz_obs = extract_Ri_second_derivative(xi, z, Ri_g)
    
    # ==========================================================================
    # STAGE 5 & 6: CLOSURE MAPPING & CURVATURE DECOMPOSITION
    # ==========================================================================
    
    # --- Track A: Businger-Dyer Closure ---
    zeta_BD = zeros(N)
    sigma_zeta_BD = zeros(N)
    kappa_zeta_BD = zeros(N)
    
    for i in 1:N
        Ri_val = Ri_g[i]
        Ri_guarded = min(Ri_val, Ri_guard)
        kappa_zeta_BD[i] = 1.0 / ((1.0 - 5.0 * Ri_guarded)^2)
        
        if Ri_val >= 0.0
            zeta_BD[i] = Ri_guarded / (1.0 - 5.0 * Ri_guarded)
            dzeta_dRi = 1.0 / ((1.0 - 5.0 * Ri_guarded)^2)
        else
            zeta_BD[i] = Ri_val
            dzeta_dRi = 1.0
        end
        
        U_z_guarded = max(abs(U_z[i]), 1e-4)
        dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z_guarded^2))
        dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z[i] / (theta_ref * U_z_guarded^3))
        sigma_zeta_BD[i] = sqrt((dzeta_dtheta_z * sigma_theta_z[i])^2 + (dzeta_dU_z * sigma_U_z[i])^2)
    end
    
    # Extract derivatives of the BD stability coordinate
    zeta_z_BD, zeta_zz_BD = extract_zeta_derivatives(xi, z, zeta_BD)
    
    # Compute BD analytical constitutive curvature components:
    # Ri_zeta = 1 / (1 + 5*zeta)^2, Ri_zeta_zeta = -10 / (1 + 5*zeta)^3 for zeta >= 0
    C_const_BD = zeros(N)
    C_coord_BD = zeros(N)
    for i in 1:N
        z_val = zeta_BD[i]
        if z_val >= 0.0
            Ri_zeta = 1.0 / ((1.0 + 5.0 * z_val)^2)
            Ri_zeta_zeta = -10.0 / ((1.0 + 5.0 * z_val)^3)
        else
            Ri_zeta = 1.0
            Ri_zeta_zeta = 0.0
        end
        C_const_BD[i] = Ri_zeta_zeta * (zeta_z_BD[i]^2)
        C_coord_BD[i] = Ri_zeta * zeta_zz_BD[i]
    end
    Ri_zz_GSPT_BD = C_const_BD .+ C_coord_BD
    delta_Ri_zz_BD = Ri_zz_obs .- Ri_zz_GSPT_BD
    
    # --- Track B: Grachev et al. (2007) Closure ---
    zeta_G07 = zeros(N)
    sigma_zeta_G07 = zeros(N)
    kappa_zeta_G07 = zeros(N)
    
    for i in 1:N
        Ri_val = Ri_g[i]
        z_sol, R_prime = invert_grachev(Ri_val)
        zeta_G07[i] = z_sol
        kappa_zeta_G07[i] = 1.0 / R_prime
        
        U_z_guarded = max(abs(U_z[i]), 1e-4)
        dzeta_dtheta_z = (1.0 / R_prime) * (g / (theta_ref * U_z_guarded^2))
        dzeta_dU_z = (1.0 / R_prime) * (-2.0 * g * theta_z[i] / (theta_ref * U_z_guarded^3))
        sigma_zeta_G07[i] = sqrt((dzeta_dtheta_z * sigma_theta_z[i])^2 + (dzeta_dU_z * sigma_U_z[i])^2)
    end
    
    # Extract derivatives of the G07 stability coordinate
    zeta_z_G07, zeta_zz_G07 = extract_zeta_derivatives(xi, z, zeta_G07)
    
    # Compute G07 analytical constitutive curvature components
    C_const_G07 = zeros(N)
    C_coord_G07 = zeros(N)
    for i in 1:N
        z_val = zeta_G07[i]
        _, Ri_zeta = R_and_prime(z_val)
        Ri_zeta_zeta = R_second_derivative_approx(z_val)
        
        C_const_G07[i] = Ri_zeta_zeta * (zeta_z_G07[i]^2)
        C_coord_G07[i] = Ri_zeta * zeta_zz_G07[i]
    end
    Ri_zz_GSPT_G07 = C_const_G07 .+ C_coord_G07
    delta_Ri_zz_G07 = Ri_zz_obs .- Ri_zz_GSPT_G07
    
    # High smoothing checks
    eta_theta = alpha_theta_irls / alpha_theta_base
    eta_U = alpha_U_irls / alpha_U_base
    is_high_smoothing = (eta_theta > 3.0 || eta_U > 3.0)
    
    return ClosureComparisonResult(
        z,
        theta_smooth, U_smooth, theta_z, U_z, Ri_g, Ri_zz_obs,
        zeta_BD, sigma_zeta_BD, kappa_zeta_BD, zeta_z_BD, zeta_zz_BD, C_const_BD, C_coord_BD, Ri_zz_GSPT_BD, delta_Ri_zz_BD,
        zeta_G07, sigma_zeta_G07, kappa_zeta_G07, zeta_z_G07, zeta_zz_G07, C_const_G07, C_coord_G07, Ri_zz_GSPT_G07, delta_Ri_zz_G07,
        is_high_smoothing, eta_theta, eta_U
    )
end

# ==============================================================================
# Numerical Demonstration & Verification Block
# ==============================================================================

# Define a 10-level logarithmically spaced SBL tower grid representing SHEBA Campaign
z_sheba = [0.5, 1.0, 1.8, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0, 14.0]

# Generate synthetic SBL profiles containing evolving stable stratification
# and a descending Low-Level Jet (LLJ) structure to trigger rich spatial curvatures
theta_true = 285.0 .+ 2.2 .* log.(z_sheba ./ 0.1) .- 0.04 .* z_sheba
# A jet structure: wind speed has a maximum near z = 8.0m
U_true = 0.5 .* log.(z_sheba ./ 0.1) .+ 0.25 .* z_sheba .- 0.015 .* (z_sheba .^ 2)

# Set realistic meteorological tower sensor noise standard deviations
delta_theta = 0.05  # Fine-wire thermocouple resolution floor (K)
delta_U = 0.02      # Sonic anemometer wind vector component resolution floor (m/s)

# Add Gaussian measurement noise for realistic verification
using Random
Random.seed!(101)
theta_obs = theta_true .+ delta_theta .* randn(length(z_sheba))
U_obs = U_true .+ delta_U .* randn(length(z_sheba))

# Run the complete SBL closure comparison pipeline
comparison = compare_sbl_closures(z_sheba, theta_obs, U_obs, delta_theta, delta_U)

# Display a high-signal scientific report of the results
println("=====================================================================================")
println("SBL PROCESSING ENGINE: STAGE 6 CLOSURE COMPARISON AUDIT")
println("=====================================================================================")
@printf("High-Smoothing Instability Flag:   %-15s\\n", comparison.is_high_smoothing ? "FLAGGED (Intermittency Suspect)" : "NORMAL")
@printf("Smoothing Scale Ratio (Theta):     %-15.2f | Smoothing Scale Ratio (Wind Speed): %.2f\\n", comparison.eta_theta, comparison.eta_U)
println("=====================================================================================")
println("STABILITY PARAMETER & COORDINATE COUPLING COMPARISON:")
println("-------------------------------------------------------------------------------------")
@printf("%-6s | %-8s | %-10s | %-10s | %-10s | %-10s | %-10s\\n", 
        "z (m)", "Ri_g", "zeta_BD", "zeta_G07", "Delta_zeta", "kappa_BD", "kappa_G07")
println("-------------------------------------------------------------------------------------")
for i in 1:length(z_sheba)
    d_zeta = comparison.zeta_BD[i] - comparison.zeta_G07[i]
    @printf("%-6.1f | %-8.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\\n",
            comparison.z[i], comparison.Ri_g[i], comparison.zeta_BD[i], comparison.zeta_G07[i], 
            d_zeta, comparison.kappa_zeta_BD[i], comparison.kappa_zeta_G07[i])
end
println("=====================================================================================")
println("SPATIAL PROFILE CURVATURE DECOMPOSITION RESIDUALS:")
println("-------------------------------------------------------------------------------------")
@printf("%-6s | %-8s | %-10s | %-10s | %-10s | %-10s | %-10s\\n",
        "z (m)", "Ri_zz_obs", "C_const_BD", "C_coord_BD", "Delta_Ri_BD", "C_const_G07", "C_coord_G07")
println("-------------------------------------------------------------------------------------")
for i in 1:length(z_sheba)
    @printf("%-6.1f | %-8.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f | %-10.4f\\n",
            comparison.z[i], comparison.Ri_zz_obs[i], 
            comparison.C_const_BD[i], comparison.C_coord_BD[i], comparison.delta_Ri_zz_BD[i],
            comparison.C_const_G07[i], comparison.C_coord_G07[i])
end
println("=====================================================================================")
