#!/usr/bin/env julia
# =============================================================================
# GSPT CURVATURE DECOMPOSITION & FOLD ILLUSION SIMULATION UTILITY
# Developed for Generalized Similarity Profile Theory (GSPT) Campaign Audits
# =============================================================================
# This script implements the mathematical and numerical foundations of GSPT:
# 1. Non-uniform finite-difference stencils (D1 and D2) to handle CASES-99
#    irregular tower spacing (1.5 m to 55 m above ground).
# 2. Synthetic nocturnal Low-Level Jet (LLJ) wind shear & flux divergence
#    profiles to model coordinate compression (zeta_zz) near the jet core.
# 3. Richardson profile curvature partitioning:
#    Total Curvature (Ri_zz) = Intrinsic Curvature (C_const) + Coordinate Curvature (C_coord)
# 4. Pre-diagnostic smoothing (Track A Regularization) via Tikhonov filtering.
# =============================================================================

using LinearAlgebra
using Printf

# Try to load Plots dynamically for headless execution
try
    using Plots
catch e
    @warn "Plots.jl not found. The script will run and save numerical data but skip rendering."
end

# =============================================================================
# 1. NON-UNIFORM GRID FINITE-DIFFERENCE OPERATORS
# =============================================================================

"""
    stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)

Computes finite-difference weights at `z0` using arbitrary coordinates `z_stencil`
for a derivative of order `m`. Formulated by solving the local Taylor-Vandermonde system.
"""
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [Float64(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p)
    b[m + 1] = 1.0  # Select target derivative order
    return A \ b
end

"""
    build_operators(z::Vector{Float64})

Generates non-uniform first (D1) and second (D2) derivative matrices on grid `z`.
Boundary points utilize second-order asymmetric stencils to ensure uniform error convergence.
"""
function build_operators(z::Vector{Float64})
    n = length(z)
    D1 = zeros(n, n)
    D2 = zeros(n, n)

    for i in 1:n
        if i == 1
            idx = [1, 2, 3]
        elseif i == n
            idx = [n-2, n-1, n]
        else
            idx = [i-1, i, i+1]
        end

        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end

    return D1, D2
end

# =============================================================================
# 2. EMPOST STABILITY CLOSURE AND ANALYTICAL DERIVATIVES
# =============================================================================

# Gradient Richardson closure function: Ri = ζ * (1 + β_h * ζ) / (1 + β_m * ζ)^2
Ri_model(ζ, β_m, β_h) = ζ * (1.0 + β_h * ζ) / (1.0 + β_m * ζ)^2

# First derivative: R' = ∂Ri / ∂ζ
function Ri_zeta(ζ, β_m, β_h)
    return (1.0 + (2.0 * β_h - β_m) * ζ) / (1.0 + β_m * ζ)^3
end

# Second derivative: R'' = ∂²Ri / ∂ζ²
function Ri_zetazeta(ζ, β_m, β_h)
    num = 2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ
    return num / (1.0 + β_m * ζ)^4
end

# =============================================================================
# 3. PLOTTING ROUTINE
# =============================================================================

"""
    plot_gspt_decomposition(hours, z_tower, C_const, C_coord, E_error; save_path)

Renders a 3-panel side-by-side contour heatmap of the decomposed GSPT curvature components.
"""
function plot_gspt_decomposition(
    hours::Vector{Float64},
    z_tower::Vector{Float64},
    C_const::Matrix{Float64},
    C_coord::Matrix{Float64},
    E_error::Matrix{Float64};
    save_path::String = "gspt_curvature_contour.png"
)
    if !@isdefined Plots
        @warn "Plots.jl not available. Skipping PNG rendering."
        return nothing
    end

    gr() # Use high-performance GR backend

    plt_opts = (
        xlabel = "Time (Hours)",
        ylabel = "Height above ground (m)",
        yticks = (z_tower, [@sprintf("%.1f", h) for h in z_tower]),
        fill = true,
        c = :coolwarm,
        linewidth = 0.5,
        tickfont = font(9, "DejaVu Sans"),
        guidefont = font(10, "DejaVu Sans"),
        titlefont = font(11, "DejaVu Sans bold")
    )

    p1 = contourf(hours, z_tower, C_const;
                  title="Intrinsic Curvature (C_const)",
                  colorbar_title="C_const (m⁻²)",
                  plt_opts...)

    p2 = contourf(hours, z_tower, C_coord;
                  title="Coordinate Curvature (C_coord)",
                  colorbar_title="C_coord (m⁻²)",
                  plt_opts...)

    p3 = contourf(hours, z_tower, E_error;
                  title="Audit Residual (E_error)",
                  colorbar_title="E_error (m⁻²)",
                  plt_opts...)

    full_plot = plot(p1, p2, p3,
                     layout = (1, 3),
                     size = (1400, 450),
                     plot_title = "GSPT Curvature Decomposition: 12-Hour SBL Cooling Cycle (CASES-99 Geometry)",
                     plot_titlefont = font(13, "DejaVu Sans bold"),
                     margin = 5Plots.mm)

    savefig(full_plot, save_path)
    println("Successfully rendered GSPT contour panel to: $save_path")
    return full_plot
end

# =============================================================================
# 4. SIMULATION LOOP (CASES-99 Tower Footprint)
# =============================================================================

function run_gspt_simulation()
    println("=====================================================================")
    println("Initializing GSPT Curvature Decomposition (CASES-99 Geometry)")
    println("=====================================================================")

    # Irregularly spaced CASES-99 tower levels
    z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
    n_z = length(z_tower)
    D1, D2 = build_operators(z_tower)

    # 12-hour cooling cycle (50 temporal steps)
    hours = collect(range(0.0, 12.0, length=50))
    n_t = length(hours)

    # Pre-allocate output matrices
    C_const_mat = zeros(n_z, n_t)
    C_coord_mat = zeros(n_z, n_t)
    E_error_mat = zeros(n_z, n_t)

    # Site-specific grassland parameters
    β_m, β_h = 5.0, 5.0
    L0 = 20.0             # Surface Obukhov length at neutral onset (m)
    sigma_obs = 0.008     # Simulated sonic anemometer/thermocouple noise floor
    lambda_reg = 5.0      # Track A Tikhonov regularization penalty

    println("Simulating 12-hour cycle with descending Low-Level Jet...")
    
    for t_idx in 1:n_t
        t = hours[t_idx]

        # Evolving jet intensity driving vertical flux-divergence
        # Jet core descends and cooling strengthens as the night progresses
        amp = 0.1 + 0.35 * (t / 12.0)
        L_profile = [L0 * (1.0 - amp * sin(π * zi / 60.0)) for zi in z_tower]

        # Similarity stability coordinate (ζ = z / L(z))
        ζ = z_tower ./ L_profile
        ζ_z = D1 * ζ
        ζ_zz = D2 * ζ

        # Compute GSPT analytical curvature components
        ri_z = [Ri_zeta(ζ[i], β_m, β_h) for i in 1:n_z]
        ri_zz = [Ri_zetazeta(ζ[i], β_m, β_h) for i in 1:n_z]

        # Decomposition: Ri_zz = C_const + C_coord
        C_const = ri_zz .* (ζ_z .^ 2)
        C_coord = ri_z .* ζ_zz
        Ri_exact = C_const .+ C_coord

        C_const_mat[:, t_idx] = C_const
        C_coord_mat[:, t_idx] = C_coord

        # Observational environment: inject high-frequency sensor noise
        noise = sigma_obs .* sin.(t_idx .+ (1:n_z) .* 1.5)
        Ri_obs = [Ri_model(ζ[i], β_m, β_h) for i in 1:n_z] .+ noise

        # Track A Primitive/Pre-Diagnostic Tikhonov Regularization:
        # Solve the system: (I + λ_reg * D2' * D2) Ri_smooth = Ri_obs
        R = D2' * D2
        A_reg = I(n_z) + lambda_reg .* R
        Ri_smooth = A_reg \ Ri_obs

        # Evaluate the regularized spatial curvature
        M_Ri_zz = D2 * Ri_smooth

        # Compute the explicit audit residual (discretization & filtering error)
        E_error_mat[:, t_idx] = M_Ri_zz .- Ri_exact
    end

    println("Decomposition complete. Curvature matrices populated.")
    
    # Render plot if Plots.jl is available
    plot_gspt_decomposition(hours, z_tower, C_const_mat, C_coord_mat, E_error_mat;
                           save_path = "plots/gspt_curvature_contour.png")

    # Log specific audit elevations at hour 12 to demonstrate the fold illusion
    println("\n=====================================================================")
    println("GSPT NUMERICAL AUDIT AT HOUR 12 (DEEP STABLE SBL REGIME)")
    println("=====================================================================")
    println("z (m)  | Intrinsic C_const | Coordinate C_coord | Total Curvature | Error E_error")
    println("---------------------------------------------------------------------")
    for i in 1:n_z
        @printf("%5.1f  |   % 13.6f |   % 14.6f |   % 13.6f |  % 11.6f\n",
                z_tower[i], C_const_mat[i, end], C_coord_mat[i, end],
                C_const_mat[i, end] + C_coord_mat[i, end], E_error_mat[i, end])
    end
    println("=====================================================================")
    println("Physical Interpretations:")
    println("1. Near-Surface Inflection Masking (z ≈ 10m): C_coord balances C_const,")
    println("   hiding intense active dynamics under a visually straight linear profile.")
    println("2. Jet Nose Singularity (z ≈ 45m): Coordinate stretching (C_coord) dominates")
    println("   observed profile curvature, creating the visual 'Fold Illusion' knee.")
    println("=====================================================================\n")
end

# Self-run if called directly
run_gspt_simulation()
