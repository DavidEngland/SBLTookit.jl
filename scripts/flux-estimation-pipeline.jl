# ==============================================================================
# SBL Vertical Profile Processing Pipeline (Stages 1-3)
# Implements Log-Coordinate Tikhonov-Morozov Spline Smoothing 
# and Downstream Businger-Dyer Richardson-Number Inversion
# ==============================================================================

using LinearAlgebra
using Printf

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
via bisection search, targeting a chi-squared discrepancy equal to the number of levels N.
"""
function fit_with_morozov(xi::Vector{Float64}, y::Vector{Float64}, sigma::Vector{Float64})
    N = length(xi)
    low = -15.0
    high = 15.0
    tol = 1e-4
    max_iter = 100
    
    local s, s_prime, Sigma_s_prime, alpha
    for iter in 1:max_iter
        mid = (low + high) / 2.0
        alpha = 10.0^mid
        s, s_prime, Sigma_s_prime = solve_smoothing_spline(xi, y, sigma, alpha)
        
        # Weighted residual sum of squares (chi^2)
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
# Pipeline Execution Functions
# ==============================================================================

struct PipelineOutput
    z::Vector{Float64}
    theta_smooth::Vector{Float64}
    U_smooth::Vector{Float64}
    theta_z::Vector{Float64}
    U_z::Vector{Float64}
    sigma_theta_z::Vector{Float64}
    sigma_U_z::Vector{Float64}
    Ri_g::Vector{Float64}
    zeta::Vector{Float64}
    kappa_zeta::Vector{Float64}
    sigma_zeta::Vector{Float64}
    alpha_theta::Float64
    alpha_U::Float64
end

"""
    run_sbl_pipeline(z, theta_raw, U_raw, delta_theta, delta_U; g=9.81, theta_ref=290.0, Ri_guard=0.19)

Runs Stages 1 through 3 of the SBL Processing Pipeline:
1. Log-coordinate Tikhonov spline smoothing using Morozov's Discrepancy Principle.
2. Gradient Richardson number mapping with coordinate-singularity guards.
3. Inversion to Monin-Obukhov non-dimensional height zeta with analytical error propagation.
"""
function run_sbl_pipeline(
    z::Vector{Float64}, 
    theta_raw::Vector{Float64}, 
    U_raw::Vector{Float64}, 
    delta_theta::Float64, 
    delta_U::Float64; 
    g=9.81, 
    theta_ref=290.0, 
    Ri_guard=0.19
)
    N = length(z)
    
    # Stage 1: Coordinate Transformation & Spline Fitting
    # Transform physical height to dimensionless log-coordinate xi
    z_0 = 0.1 # reference surface roughness scale (treated as numerical offset)
    xi = log.(z ./ z_0)
    
    # Vectorized measurement noise floors
    sigma_theta = fill(delta_theta, N)
    sigma_U = fill(delta_U, N)
    
    # Fit potential temperature and wind speed profiles
    theta_smooth, theta_xi_prime, Sigma_theta_xi_prime, alpha_theta = fit_with_morozov(xi, theta_raw, sigma_theta)
    U_smooth, U_xi_prime, Sigma_U_xi_prime, alpha_U = fit_with_morozov(xi, U_raw, sigma_U)
    
    # Extract vertical physical derivatives: d/dz = (1/z) * d/dxi
    theta_z = theta_xi_prime ./ z
    U_z = U_xi_prime ./ z
    
    # Propagate physical gradient standard deviations
    sigma_theta_z = sqrt.(diag(Sigma_theta_xi_prime)) ./ z
    sigma_U_z = sqrt.(diag(Sigma_U_xi_prime)) ./ z
    
    # Stage 2: Gradient Richardson Number Mapping
    Ri_g = (g / theta_ref) .* theta_z ./ (U_z .^ 2)
    
    # Apply Richardson-number guard near Businger-Dyer asymptotic limit (Ri_g -> 0.2)
    Ri_guarded = min.(Ri_g, Ri_guard)
    
    # Compute ill-conditioning metric (kappa_zeta)
    kappa_zeta = 1.0 ./ ((1.0 .- 5.0 .* Ri_guarded) .^ 2)
    
    # Stage 3: Monin-Obukhov Similarity Inversion & Analytical Jacobian
    zeta = zeros(N)
    sigma_zeta = zeros(N)
    
    for i in 1:N
        # Businger-Dyer Inversion and Analytical derivative dzeta/dRi_g
        dzeta_dRi = 0.0
        if Ri_g[i] >= 0.0
            zeta[i] = Ri_guarded[i] / (1.0 - 5.0 * Ri_guarded[i])
            dzeta_dRi = 1.0 / ((1.0 - 5.0 * Ri_guarded[i])^2)
        else
            zeta[i] = Ri_g[i]
            dzeta_dRi = 1.0
        end
        
        # Chain rule sensitivities: dzeta/dtheta_z and dzeta/dU_z
        dzeta_dtheta_z = dzeta_dRi * (g / (theta_ref * U_z[i]^2))
        dzeta_dU_z = dzeta_dRi * (-2.0 * g * theta_z[i] / (theta_ref * U_z[i]^3))
        
        # Analytical first-order uncertainty propagation
        sigma_zeta[i] = sqrt(
            (dzeta_dtheta_z * sigma_theta_z[i])^2 + 
            (dzeta_dU_z * sigma_U_z[i])^2
        )
    end
    
    return PipelineOutput(
        z, theta_smooth, U_smooth, theta_z, U_z, 
        sigma_theta_z, sigma_U_z, Ri_g, zeta, kappa_zeta, 
        sigma_zeta, alpha_theta, alpha_U
    )
end

# ==============================================================================
# Numerical Demonstration & Verification Block
# ==============================================================================

# Define a 10-level logarithmically spaced meteorological tower grid (e.g., SHEBA-like)
z = [0.5, 1.0, 1.8, 3.0, 4.5, 6.0, 8.0, 10.0, 12.0, 14.0]

# True noise-free physical profiles
theta_true = 285.0 .+ 2.0 .* log.(z ./ 0.1) .- 0.05 .* z
U_true = 0.4 .* log.(z ./ 0.1) .+ 0.1 .* z

# Realistic instrument noise floors
delta_theta = 0.05 # standard fine-wire thermocouple precision (K)
delta_U = 0.02     # sonic anemometer resolution floor (m/s)

# Generate synthetic observations with Gaussian sensor noise
using Random
Random.seed!(42) # Ensure exact reproducibility of noisy data
theta_obs = theta_true .+ delta_theta .* randn(length(z))
U_obs = U_true .+ delta_U .* randn(length(z))

# Run the complete pipeline
res = run_sbl_pipeline(z, theta_obs, U_obs, delta_theta, delta_U)

# Display results
println("=====================================================================================")
println("SBL PROCESSING ENGINE: STAGES 1 TO 3 VERIFICATION")
println("Morozov Parameter Theta (alpha): ", round(res.alpha_theta, sigdigits=4))
println("Morozov Parameter Wind (alpha):  ", round(res.alpha_U, sigdigits=4))
println("=====================================================================================")
@printf("%-6s | %-8s | %-8s | %-8s | %-8s | %-8s | %-8s\n", 
        "z (m)", "Theta_z", "U_z", "Ri_g", "kappa_z", "zeta", "sigma_zeta")
println("-------------------------------------------------------------------------------------")
for i in 1:length(z)
    @printf("%-6.1f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f | %-8.4f\n",
            res.z[i], res.theta_z[i], res.U_z[i], res.Ri_g[i], res.kappa_zeta[i], res.zeta[i], res.sigma_zeta[i])
end
println("=====================================================================================")
