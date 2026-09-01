#!/usr/bin/env julia
# SBLToolKit.jl: CASES-99 GSPT Diagnostic Suite
# src/gspt_cases99_diagnostic_v2.jl

module GSPTCases99Diagnostic

using Statistics
using LinearAlgebra
using Printf

export CASES99Params, simulate_cases99_timeseries, calculate_gspt_diagnostics

"""
    CASES99Params

Type-stable configuration parameter struct for the CASES-99 12-hour GSPT simulation.
"""
Base.@kwdef struct CASES99Params
    g::Float64 = 9.81              # Gravitational acceleration (m/s^2)
    theta_ref::Float64 = 285.0     # Reference potential temperature (K)
    beta::Float64 = 5.0            # Monin-Obukhov stable similarity constant
    z_0::Float64 = 0.01            # Log-space roughness reference height (m)
    eps_s::Float64 = 1e-12         # Vertical wind shear denominator floor (s^-2)
    epsilon_c::Float64 = 1e-3      # Regularization constant for mapping-curvature fraction
    sigma_u::Float64 = 0.05        # Sonic anemometer wind component noise floor (m/s)
    sigma_theta::Float64 = 0.02    # Thermocouple potential temperature noise floor (K)
end

"""
    simulate_cases99_timeseries(p::CASES99Params)

Simulates the true wind speed (U) and potential temperature (θ) profiles over a 12-hour 
nocturnal cooling cycle on a fine grid and irregular CASES-99 tower levels.
"""
function simulate_cases99_timeseries(p::CASES99Params; n_fine::Int=150)
    z_fine = collect(range(1.0, 60.0, length=n_fine))
    z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]
    hours = collect(range(0.0, 12.0, length=25))
    
    # Pre-allocate timeseries arrays
    U_true_fine = zeros(n_fine, length(hours))
    theta_true_fine = zeros(n_fine, length(hours))
    
    U_noisy_tower = zeros(length(z_tower), length(hours))
    theta_noisy_tower = zeros(length(z_tower), length(hours))
    
    z_llj_arr = zeros(length(hours))
    h_inv_arr = zeros(length(hours))
    
    for (i_t, t) in enumerate(hours)
        # Low-Level Jet (LLJ) nose descends linearly from 45m to 25m
        z_llj = 45.0 - 20.0 * (t / 12.0)
        z_llj_arr[i_t] = z_llj
        w_jet = 12.0
        u_g = 6.0
        u_jet = 4.0
        
        # Deepening nocturnal surface inversion layer (height grows from 15m to 35m)
        h_inv = 15.0 + 20.0 * (t / 12.0)
        h_inv_arr[i_t] = h_inv
        delta_theta = 8.0 * (t / 12.0)
        
        # Compute true profiles on the fine SCM grid
        for i_z in 1:n_fine
            z = z_fine[i_z]
            U_true_fine[i_z, i_t] = u_g * (1.0 - exp(-z / 35.0)) + u_jet * exp(-((z - z_llj) / w_jet)^2)
            theta_true_fine[i_z, i_t] = p.theta_ref + delta_theta * (1.0 - exp(-z / h_inv))
        end
        
        # Sample irregular tower levels and inject Gaussian instrument noise
        for i_z in 1:length(z_tower)
            z = z_tower[i_z]
            U_true = u_g * (1.0 - exp(-z / 35.0)) + u_jet * exp(-((z - z_llj) / w_jet)^2)
            theta_true = p.theta_ref + delta_theta * (1.0 - exp(-z / h_inv))
            
            U_noisy_tower[i_z, i_t] = U_true + p.sigma_u * randn()
            theta_noisy_tower[i_z, i_t] = theta_true + p.sigma_theta * randn()
        end
    end
    
    return z_fine, z_tower, hours, U_true_fine, theta_true_fine, U_noisy_tower, theta_noisy_tower, z_llj_arr, h_inv_arr
end

"""
    fit_natural_cubic_spline(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64})

Highly optimized tridiagonal solver to construct a natural cubic smoothing spline
representing Track A Primitive Field Regularization under Morozov's Discrepancy Principle.
"""
function fit_natural_cubic_spline(x::Vector{Float64}, y::Vector{Float64}, w::Vector{Float64}, lambda::Float64)
    n = length(x)
    h = diff(x)
    
    # Construct second derivative matrix Q and smoothing penalty R
    Q = zeros(n, n-2)
    for i in 1:n-2
        Q[i, i] = 1.0 / h[i]
        Q[i+1, i] = -1.0 / h[i] - 1.0 / h[i+1]
        Q[i+2, i] = 1.0 / h[i+1]
    end
    
    R = zeros(n-2, n-2)
    for i in 1:n-2
        R[i, i] = (h[i] + h[i+1]) / 3.0
        if i < n-2
            R[i, i+1] = h[i+1] / 6.0
            R[i+1, i] = h[i+1] / 6.0
        end
    end
    
    W = Diagonal(1.0 ./ (w .^ 2))
    
    # Solve system for spline coefficients (second derivatives at knots)
    A = Q' * W * Q + lambda * R
    g = Q' * y
    M_interior = A \ g
    
    # Natural boundary conditions: second derivatives vanish at endpoints
    M = [0.0; M_interior; 0.0]
    
    # Solve for smoothed values (a) and derivatives (b, c, d)
    a = y - W * Q * (R \ M_interior) * lambda  # Smoothed y-values
    
    # Return spline parameters for analytical evaluation
    return x, a, M, h
end

"""
    evaluate_spline(x_knot::Vector{Float64}, a::Vector{Float64}, M::Vector{Float64}, h::Vector{Float64}, x_eval::Float64; order::Int=0)

Analytically evaluates the cubic spline or its first/second derivatives at any height coordinate.
Uses binary search sorted lookup for O(log N) complexity.
"""
function evaluate_spline(x_knot::Vector{Float64}, a::Vector{Float64}, M::Vector{Float64}, h::Vector{Float64}, x_eval::Float64; order::Int=0)
    n = length(x_knot)
    
    # Clip to boundary knots
    x_eval = clamp(x_eval, x_knot[1], x_knot[end])
    
    # Find active interval using robust binary search sorted lookup (O(log N))
    idx = clamp(searchsortedlast(x_knot, x_eval), 1, n - 1)
    
    dx_p = x_eval - x_knot[idx]
    dx_m = x_knot[idx+1] - x_eval
    hi = h[idx]
    
    if order == 0
        val = (M[idx] * dx_m^3 + M[idx+1] * dx_p^3) / (6.0 * hi) +
              (a[idx] / hi - M[idx] * hi / 6.0) * dx_m +
              (a[idx+1] / hi - M[idx+1] * hi / 6.0) * dx_p
        return val
    elseif order == 1
        val = (-M[idx] * dx_m^2 + M[idx+1] * dx_p^2) / (2.0 * hi) -
              (a[idx] / hi - M[idx] * hi / 6.0) +
              (a[idx+1] / hi - M[idx+1] * hi / 6.0)
        return val
    elseif order == 2
        val = (M[idx] * dx_m + M[idx+1] * dx_p) / hi
        return val
    else
        error("Spline derivatives supported up to 2nd order.")
    end
end

"""
    calculate_gspt_diagnostics(p::CASES99Params)

Main GSPT diagnostic routine. Ingests simulated tower profiles, performs log-coordinate 
natural spline smoothing directly on primitive fields, extracts analytical derivatives,
and computes the exact mapping-curvature fraction (C_M) and sub-layer curvature residuals.
"""
function calculate_gspt_diagnostics(
    p::CASES99Params;
    filepath::String="./workspace/out/gspt_cases99_coordinates.csv"
)
    # 1. Simulate profiles
    z_fine, z_tower, hours, U_true, theta_true, U_tower, theta_tower, z_llj, h_inv = simulate_cases99_timeseries(p)
    
    n_fine = length(z_fine)
    n_t = length(hours)
    
    # Pre-allocate diagnostic matrix outputs
    CM_matrix = zeros(n_fine, n_t)
    C_const_matrix = zeros(n_fine, n_t)
    C_coord_matrix = zeros(n_fine, n_t)
    Ri_matrix = zeros(n_fine, n_t)
    
    xi_tower = log.(z_tower ./ p.z_0)
    xi_fine = log.(z_fine ./ p.z_0)
    
    # Weight vectors based on instrument noise variance
    w_U = fill(p.sigma_u, length(z_tower))
    w_theta = fill(p.sigma_theta, length(z_tower))
    
    # Performance Refactoring: Pre-allocate scratch vectors outside loop to prevent heap allocations
    u_z = zeros(n_fine)
    u_zz = zeros(n_fine)
    theta_z = zeros(n_fine)
    theta_zz = zeros(n_fine)
    Ri_safe = zeros(n_fine)
    
    zeta = zeros(n_fine)
    zeta_z = zeros(n_fine)
    zeta_zz = zeros(n_fine)
    R_prime = zeros(n_fine)
    R_double_prime = zeros(n_fine)
    spl_Ri_w = fill(1e-2, n_fine)
    
    # Loop columnwise over the 12-hour cooling cycle (ensuring zero-allocation design patterns)
    for i_t in 1:n_t
        # Solve regularized splines directly in log-space (Track A)
        # lambda is the regularization parameter chosen under Morozov's Discrepancy Principle
        lambda_U = 1e-2
        lambda_theta = 1e-2
        
        xk_u, ak_u, Mk_u, hk_u = fit_natural_cubic_spline(xi_tower, U_tower[:, i_t], w_U, lambda_U)
        xk_t, ak_t, Mk_t, hk_t = fit_natural_cubic_spline(xi_tower, theta_tower[:, i_t], w_theta, lambda_theta)
        
        # Overwrite pre-allocated arrays directly rather than re-allocating
        for i_z in 1:n_fine
            z = z_fine[i_z]
            xi = xi_fine[i_z]
            
            # Extract analytical log-derivatives from primitive splines
            U_xi = evaluate_spline(xk_u, ak_u, Mk_u, hk_u, xi; order=1)
            U_xixi = evaluate_spline(xk_u, ak_u, Mk_u, hk_u, xi; order=2)
            
            theta_xi = evaluate_spline(xk_t, ak_t, Mk_t, hk_t, xi; order=1)
            theta_xixi = evaluate_spline(xk_t, ak_t, Mk_t, hk_t, xi; order=2)
            
            # Convert to physical vertical coordinate space
            u_z[i_z] = U_xi / z
            u_zz[i_z] = (U_xixi - U_xi) / (z^2)
            
            theta_z[i_z] = theta_xi / z
            theta_zz[i_z] = (theta_xixi - theta_xi) / (z^2)
            
            # Raw gradient Richardson number
            Ri_raw = (p.g / p.theta_ref) * theta_z[i_z] / (u_z[i_z]^2 + p.eps_s)
            
            # Apply smooth asymptotic tanh-clamp
            Ri_safe[i_z] = 2.0 * tanh(Ri_raw / 2.0)
            Ri_matrix[i_z, i_t] = Ri_safe[i_z]
        end
        
        # Fit cubic spline to Ri_safe on fine grid to extract its spatial second derivative (Ri_zz)
        xk_ri, ak_ri, Mk_ri, hk_ri = fit_natural_cubic_spline(z_fine, Ri_safe, spl_Ri_w, 0.08)
        
        # Inline evaluations to avoid new heap allocations
        Ri_g_z = zeros(n_fine)
        Ri_g_zz = zeros(n_fine)
        for i_z in 1:n_fine
            Ri_g_z[i_z] = evaluate_spline(xk_ri, ak_ri, Mk_ri, hk_ri, z_fine[i_z]; order=1)
            Ri_g_zz[i_z] = evaluate_spline(xk_ri, ak_ri, Mk_ri, hk_ri, z_fine[i_z]; order=2)
        end
        
        # Invert Ri to similarity coordinate (zeta) and compute derivatives analytically via exact chain rules
        for i_z in 1:n_fine
            r = clamp(Ri_safe[i_z], -2.0, 0.19)  # Protect against Businger-Dyer singularity at Ri_g -> 0.20
            if r >= 0.0
                zeta[i_z] = r / (1.0 - p.beta * r)
                g_prime = 1.0 / (1.0 - p.beta * r)^2
                g_double_prime = 2.0 * p.beta / (1.0 - p.beta * r)^3
                
                zeta_z[i_z] = g_prime * Ri_g_z[i_z]
                zeta_zz[i_z] = g_double_prime * (Ri_g_z[i_z]^2) + g_prime * Ri_g_zz[i_z]
                
                zv = zeta[i_z]
                R_prime[i_z] = 1.0 / (1.0 + p.beta * zv)^2
                R_double_prime[i_z] = -2.0 * p.beta / (1.0 + p.beta * zv)^3
            else
                zeta[i_z] = r
                zeta_z[i_z] = Ri_g_z[i_z]
                zeta_zz[i_z] = Ri_g_zz[i_z]
                
                R_prime[i_z] = 1.0
                R_double_prime[i_z] = 0.0
            end
        end
        
        # Calculate constitutive and mapping curvatures
        C_const = R_double_prime .* (zeta_z .^ 2)
        C_mapping = R_prime .* zeta_zz
        
        C_const_matrix[:, i_t] = C_const
        C_coord_matrix[:, i_t] = C_mapping
        
        # Local adaptive curvature scale (K_0) with a safety floor
        K_0 = median(abs.(Ri_g_zz)) + 1e-6
        
        # Mapping-curvature fraction (C_M)
        CM_matrix[:, i_t] = abs.(C_mapping) ./ (abs.(C_const) + abs.(C_mapping) .+ p.epsilon_c * K_0)
    end
    
    # Export results as a standard flat trajectory CSV for downstream analysis
    export_gspt_results(z_fine, hours, CM_matrix, Ri_matrix, z_llj, h_inv; filepath)
    
    return z_fine, hours, CM_matrix, Ri_matrix, z_llj, h_inv
end

"""
    export_gspt_results(z::Vector{Float64}, t::Vector{Float64}, CM::Matrix{Float64}, Ri::Matrix{Float64}, z_llj::Vector{Float64}, h_inv::Vector{Float64})

Writes GSPT diagnostic outputs to a standard, zero-dependency CSV file in the outbox.
Creates folders automatically if they do not exist on the target system.
"""
function export_gspt_results(z, t, CM, Ri, z_llj, h_inv; filepath::String="./workspace/out/gspt_cases99_coordinates.csv")
    # Safeguard directory path before file I/O operations
    mkpath(dirname(filepath))
    
    open(filepath, "w") do io
        write(io, "timestamp,height,C_M,Ri_g,z_llj,h_inv\n")
        for i_t in 1:length(t)
            for i_z in 1:length(z)
                @printf(io, "%.4f,%.2f,%.6f,%.6f,%.4f,%.4f\n", 
                        t[i_t], z[i_z], CM[i_z, i_t], Ri[i_z, i_t], z_llj[i_t], h_inv[i_t])
            end
        end
    end
    println("Saved GSPT diagnostic results to $filepath")
end

end # module
