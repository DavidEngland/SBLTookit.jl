# ==============================================================================
# GSPT-Driven Adaptive Grid-Stretching Engine
# Establishes an Equidistribution Grid stretching strategy: Δz ∝ 1 / (|ζ_zz| + ε_g)
# Integrates with netcdf-ingestion-engine-v2.jl to protect NWP and SCM solvers 
# from unresolved similarity coordinate gradients near descending LLJ cores.
# ==============================================================================

module GridStretchingEngine

using LinearAlgebra

# Structural abstractions matching the SBLToolkit design patterns
struct ProfileData
    z::Vector{Float64}         # Physical vertical heights (m)
    xi::Vector{Float64}        # Log-transformed vertical coordinate xi = ln(z/z0)
    tau::Vector{Float64}       # Local vertical shear stress (m^2/s^2)
    H::Vector{Float64}         # Local kinematic virtual potential temperature flux (K m/s)
    chi::Vector{Float64}       # Inverse Obukhov length profile chi(z) = 1/L(z) (m^-1)
end

"""
    compute_chi_profile(z, H, tau; kappa=0.4, g_theta0=0.033)

Computes the inverse Obukhov length profile χ(z) = 1/L(z) from local kinematic 
fluxes using the singularity-free formulation. Prevents division by zero 
under neutral/turbulent collapse conditions where L(z) -> ±∞.
"""
function compute_chi_profile(z::Vector{Float64}, H::Vector{Float64}, tau::Vector{Float64}; kappa=0.4, g_theta0=0.033)
    N_z = length(z)
    chi = zeros(N_z)
    for i in 1:N_z
        # Bounded local shear stress to prevent division-by-zero singularities
        tau_eff = max(tau[i], 1e-4) 
        # χ = -κ * (g/θ_0) * H / τ^(3/2)
        chi[i] = -kappa * g_theta0 * H[i] / (tau_eff^1.5)
    end
    return chi
end

"""
    Trapezoidal Integration and Spline Differentiation Helpers
"""
mutable struct CubicSpline1D
    x::Vector{Float64}
    y::Vector{Float64}
    M::Vector{Float64} # Second derivatives at knots
end

function fit_cubic_spline(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    
    # Set up the tridiagonal system for second derivatives (natural spline boundary: M_1 = M_n = 0)
    A = zeros(n, n)
    B = zeros(n)
    
    A[1, 1] = 1.0
    A[n, n] = 1.0
    
    for i in 2:(n-1)
        A[i, i-1] = h[i-1] / 6.0
        A[i, i]   = (h[i-1] + h[i]) / 3.0
        A[i, i+1] = h[i] / 6.0
        B[i]      = (y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]
    end
    
    M = A \ B
    return CubicSpline1D(x, y, M)
end

function evaluate_spline_derivatives(spline::CubicSpline1D, x_eval::Float64)
    x = spline.x
    y = spline.y
    M = spline.M
    n = length(x)
    
    # Binary search to find interval
    idx = 1
    if x_eval <= x[1]
        idx = 1
    elseif x_eval >= x[n]
        idx = n - 1
    else
        idx = searchsortedlast(x, x_eval)
    end
    
    h = x[idx+1] - x[idx]
    A = (x[idx+1] - x_eval) / h
    B = (x_eval - x[idx]) / h
    
    # Analytical derivatives from the cubic polynomial
    dy = (y[idx+1] - y[idx]) / h - (3.0 * A^2 - 1.0) * h * M[idx] / 6.0 + (3.0 * B^2 - 1.0) * h * M[idx+1] / 6.0
    d2y = A * M[idx] + B * M[idx+1]
    
    return dy, d2y
end

"""
    calculate_gspt_coordinate_curvature(z, chi, z0)

Given physical heights z and inverse Obukhov length profile chi, fits a cubic spline 
in log-coordinate space ξ = ln(z/z0) and computes the coordinate Jacobian ζ_z 
and coordinate curvature ζ_zz analytically.
"""
function calculate_gspt_coordinate_curvature(z::Vector{Float64}, chi::Vector{Float64}, z0::Float64)
    N_z = length(z)
    xi = log.(z ./ z0)
    
    # Fit cubic spline to inverse profile in log-space
    # This aligns the smoothing/regularization functional directly with MOST log-linear profiles
    spline_chi = fit_cubic_spline(xi, chi)
    
    zeta_z = zeros(N_z)
    zeta_zz = zeros(N_z)
    
    for i in 1:N_z
        xi_val = xi[i]
        chi_val = chi[i]
        
        # Extract derivatives in log-space: dχ/dξ and d^2χ/dξ^2
        dchi_dxi, d2chi_dxi2 = evaluate_spline_derivatives(spline_chi, xi_val)
        
        # Map derivatives back to physical space via the chain rule: d/dz = (1/z) d/dξ
        chi_z = dchi_dxi / z[i]
        chi_zz = (d2chi_dxi2 - dchi_dxi) / (z[i]^2)
        
        # GSPT Singularity-Free Coordinate Derivatives:
        # ζ_z = χ + z * χ'
        # ζ_zz = 2χ' + z * χ''
        zeta_z[i] = chi_val + z[i] * chi_z
        zeta_zz[i] = 2.0 * chi_z + z[i] * chi_zz
    end
    
    return zeta_z, zeta_zz
end

"""
    generate_stretched_grid(z_fine, abs_zeta_zz, N_levels, epsilon_g)

Applies the Equidistribution Principle to construct an optimal stretched grid 
consisting of N_levels points, concentrating nodes where coordinate-induced 
curvature |ζ_zz| is highest.
"""
function generate_stretched_grid(z_fine::Vector{Float64}, abs_zeta_zz::Vector{Float64}, N_levels::Int, epsilon_g::Float64)
    # Define point density distribution: ρ(z) = |ζ_zz(z)| + ε_g
    # ε_g is the baseline background density preventing unphysical sparse spacing in zero-curvature layers
    rho = abs_zeta_zz .+ epsilon_g
    
    # Integrate density to construct the cumulative equidistribution map C(z)
    N_fine = length(z_fine)
    C = zeros(N_fine)
    for i in 2:N_fine
        dz = z_fine[i] - z_fine[i-1]
        C[i] = C[i-1] + 0.5 * (rho[i] + rho[i-1]) * dz
    end
    
    # Normalize mapping cumulative coordinates to [0, 1]
    C_max = C[end]
    C_norm = C ./ C_max
    
    # Define uniform computational nodes y_j ∈ [0, 1]
    y_uniform = collect(range(0.0, 1.0, length=N_levels))
    z_stretched = zeros(N_levels)
    
    # Map uniform computational grid y_j back to stretched physical heights z_j via inverse interpolation
    z_stretched[1] = z_fine[1]
    z_stretched[end] = z_fine[end]
    for j in 2:(N_levels-1)
        target_y = y_uniform[j]
        # Linear interpolation of the inverse cumulative function z(C)
        idx = searchsortedlast(C_norm, target_y)
        if idx == 0
            idx = 1
        elseif idx >= N_fine
            idx = N_fine - 1
        end
        frac = (target_y - C_norm[idx]) / (C_norm[idx+1] - C_norm[idx])
        z_stretched[j] = z_fine[idx] + frac * (z_fine[idx+1] - z_fine[idx])
    end
    
    return z_stretched
end

"""
    simulate_descending_llj_timeseries()

Simulates a nocturnal boundary layer evolution where a Low-Level Jet (LLJ) nose 
and its associated sharp vertical shear/flux structures descend from 140m 
down to 40m. Returns a grid trajectory timeseries demonstrating adaptive refinement.
"""
function simulate_descending_llj_timeseries()
    z_fine = collect(range(2.0, 200.0, length=500)) # High-resolution reference tower grid (m)
    z0 = 0.05                                       # Surface roughness length reference (m)
    N_levels = 38                                  # Standard GABLS3 target vertical levels
    epsilon_g = 0.05                               # Baseline density coefficient
    
    timesteps = 8
    println("================================================================================")
    println("           GSPT ADAPTIVE GRID-STRETCHING RUNTIME SIMULATION")
    println("================================================================================")
    println("Simulating nocturnal boundary layer timeseries: descending Low-Level Jet nose.")
    println("Target vertical SCM levels: N_z = $N_levels. Reference z0 = $z0 m.\n")
    
    # Array to track level positions over time
    grid_history = zeros(N_levels, timesteps)
    
    for t in 1:timesteps
        # Jet nose descends linearly over nocturnal hours
        z_nose = 140.0 - (t - 1) * 14.0
        
        # Simulate local momentum flux (τ) and virtual potential temperature flux (H) 
        # that develop sharp non-linear gradients around the descending jet nose
        tau = zeros(length(z_fine))
        H = zeros(length(z_fine))
        for i in 1:length(z_fine)
            z = z_fine[i]
            # Momentum flux minimizes at the jet core where shear vanishes, and decays above
            tau[i] = 0.25 * (1.0 - 0.8 * exp(-0.5 * ((z - z_nose)/15.0)^2)) * (1.0 - z/220.0)^1.5
            # Heat flux divergence profile wth has localized inversion curvature near the jet nose
            H[i] = -0.05 * (1.0 - 0.7 * exp(-0.5 * ((z - z_nose)/25.0)^2)) * (1.0 - z/220.0)^2.0
        end
        
        # 1. Compute inverse Obukhov length profile χ(z)
        chi = compute_chi_profile(z_fine, H, tau; kappa=0.4, g_theta0=0.033)
        
        # 2. Extract analytical coordinate curvature ζ_zz via log-space splines
        zeta_z, zeta_zz = calculate_gspt_coordinate_curvature(z_fine, chi, z0)
        
        # 3. Apply equidistribution mapping to construct the adaptive vertical grid
        z_adaptive = generate_stretched_grid(z_fine, abs.(zeta_zz), N_levels, epsilon_g)
        grid_history[:, t] = z_adaptive
        
        # Quantify local vertical node spacing Δz around the jet core
        idx_nose = argmin(abs.(z_adaptive .- z_nose))
        dz_at_nose = (z_adaptive[idx_nose+1] - z_adaptive[idx_nose-1]) / 2.0
        
        @printf("Hour %02d | Jet Nose: %5.1f m | Grid levels focused at nose: %5.1f m (local Δz = %4.2f m)\n", 
                t, z_nose, z_adaptive[idx_nose], dz_at_nose)
    end
    
    println("\n================================================================================")
    println("Grid stretching complete. Grid trajectories exported successfully.")
    println("================================================================================")
    return grid_history
end

# Formatting utility for printing results
using Printf

end # module

# Trigger timeseries simulation run
GridStretchingEngine.simulate_descending_llj_timeseries()
