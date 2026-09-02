using LinearAlgebra
using Printf

# ==============================================================================
# GSPT Spatial Operators & Curvature Audit Engine
# ==============================================================================

"""
    stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)

Computes the finite-difference stencil weights at a target coordinate `z0` using
an arbitrary grid set `z_stencil` for a derivative of order `m`.
This is formulated by solving a Vandermonde-like system derived from local Taylor expansion.
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

Generates non-uniform derivative matrices D1 (first derivative) and D2 (second derivative)
on a vertical coordinate vector `z`. Boundary points utilize second-order accurate asymmetric
one-sided stencils to ensure uniform error convergence across all levels.
"""
function build_operators(z::Vector{Float64})
    n = length(z)
    D1 = zeros(n, n)
    D2 = zeros(n, n)

    for i in 1:n
        if i == 1
            # Lower boundary: Asymmetric one-sided 3-point stencil
            idx = [1, 2, 3]
        elseif i == n
            # Upper boundary: Asymmetric one-sided 3-point stencil
            idx = [n-2, n-1, n]
        else
            # Interior levels: Non-uniform centered 3-point stencil
            idx = [i-1, i, i+1]
        end

        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end

    return D1, D2
end

# Constitutive stability derivatives (General Case: β_h != β_m)
Ri_model(ζ, β_m, β_h) = ζ * (1.0 + β_h * ζ) / (1.0 + β_m * ζ)^2

function Ri_zeta(ζ, β_m, β_h)
    return (1.0 + (2.0 * β_h - β_m) * ζ) / (1.0 + β_m * ζ)^3
end

function Ri_zetazeta(ζ, β_m, β_h)
    num = 2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ
    return num / (1.0 + β_m * ζ)^4
end

"""
    audit_gspt_curvature(z, Ri_obs, L_profile, β_m, β_h)

Performs the three-layer GSPT curvature audit along a non-uniform grid.
Separates intrinsic thermodynamic curvature (C_const) from coordinate stretching (C_coord)
and isolates the discrete noncommutation/discretization residual (E_error) explicitly.
"""
function audit_gspt_curvature(
    z::Vector{Float64},
    Ri_obs::Vector{Float64},
    L_profile::Vector{Float64},
    β_m::Float64,
    β_h::Float64
)
    n = length(z)
    D1, D2 = build_operators(z)

    # Calculate similarity coordinates: ζ = z/L(z)
    ζ_true = z ./ L_profile

    # Evaluate coordinate spatial derivatives on the SCM grid
    ζ_z = D1 * ζ_true
    ζ_zz = D2 * ζ_true

    # Evaluate analytical constitutive derivatives
    Ri_z_val = [Ri_zeta(ζ_true[i], β_m, β_h) for i in 1:n]
    Ri_zz_val = [Ri_zetazeta(ζ_true[i], β_m, β_h) for i in 1:n]

    # Decompose Curvature Layers
    C_const = Ri_zz_val .* (ζ_z .^ 2)        # Intrinsic stability geometry
    C_coord = Ri_z_val .* ζ_zz              # Flux-coordinate geometry
    Ri_exact = C_const .+ C_coord           # Total continuous analytical curvature

    # Calculate observed physical curvature of the profile
    # (subject to discrete truncation / coarse-graining bias)
    Ri_zz_obs = D2 * Ri_obs

    # Isolate discrete audit residual (truncation/noncommutation error)
    E_error = Ri_zz_obs .- Ri_exact

    return C_const, C_coord, Ri_exact, Ri_zz_obs, E_error
end