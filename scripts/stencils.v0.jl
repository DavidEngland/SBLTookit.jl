using LinearAlgebra, Printf, Random

# --- 1. Non-Uniform Finite-Difference Stencil Weights ---
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p); b[m + 1] = 1.0
    return A \ b
end

function build_operators(z::Vector{Float64})
    n = length(z)
    D1, D2 = zeros(n, n), zeros(n, n)

    # Interior 3-point stencils
    for i in 2:n-1
        idx = [i-1, i, i+1]
        D1[i, idx] .= stencil_weights(z[idx], z[i], 1)
        D2[i, idx] .= stencil_weights(z[idx], z[i], 2)
    end

    # One-sided 3-point boundary stencils
    D1[1, 1:3] .= stencil_weights(z[1:3], z[1], 1)
    D2[1, 1:3] .= stencil_weights(z[1:3], z[1], 2)
    D1[n, n-2:n] .= stencil_weights(z[n-2:n], z[n], 1)
    D2[n, n-2:n] .= stencil_weights(z[n-2:n], z[n], 2)

    return D1, D2
end

# --- 2. Asymmetric Stability Functions & Derivatives ---
# General case: Ri(ζ) = ζ(1 + β_h*ζ) / (1 + β_m*ζ)²
Ri_model(ζ, β_m, β_h) = ζ * (1.0 + β_h * ζ) / (1.0 + β_m * ζ)^2

function Ri_zeta(ζ, β_m, β_h)
    return (1.0 + (2.0 * β_h - β_m) * ζ) / (1.0 + β_m * ζ)^3
end

function Ri_zetazeta(ζ, β_m, β_h)
    num = 2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ
    return num / (1.0 + β_m * ζ)^4
end

# --- 3. CASES-99 Execution Setup ---
z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
n = length(z_tower)
D1, D2 = build_operators(z_tower)

# Model parameters for asymmetric nocturnal SBL
β_m, β_h = 1.0, 3.0
L0 = 20.0  # Base Obukhov length (m)

# Nocturnal LLJ Flux Divergence Profile: L(z) collapses near LLJ nose (z ≈ 35 m)
L_profile(z) = L0 * (1.0 - 0.45 * sin(π * z / 60.0))

# Evaluate continuous similarity coordinates and coordinate derivatives
ζ_true = [z / L_profile(z) for z in z_tower]
ζ_z    = D1 * ζ_true
ζ_zz   = D2 * ζ_true

# Analytical Curvature Components
Ri_z_val  = [Ri_zeta(ζ_true[i], β_m, β_h) for i in 1:n]
Ri_zz_val = [Ri_zetazeta(ζ_true[i], β_m, β_h) for i in 1:n]

C_const = Ri_zz_val .* (ζ_z .^ 2)        # Constitutive similarity geometry
C_coord = Ri_z_val .* ζ_zz              # Flux-coordinate geometry
Ri_exact = C_const .+ C_coord           # Total continuous analytical curvature

# --- 4. Tikhonov Regularized Observation Operator (M_Δz) ---
Random.seed!(42)
σ_obs = 0.008
Ri_obs = [Ri_model(ζ_true[i], β_m, β_h) for i in 1:n] .+ σ_obs .* randn(n)

# Tikhonov smoother: (I + λ D₂ᵀ D₂) Ri_smooth = Ri_obs
R = D2' * D2
λ = 5.0  # Morozov/GCV regularization parameter
Ri_smooth = (I(n) + λ .* R) \ Ri_obs

# Observed Physical Curvature
M_Ri_zz = D2 * Ri_smooth

# Estimation / Discretization Error
E_error = M_Ri_zz .- Ri_exact

# --- 5. Display GSPT Curvature Strata Breakdown ---
println("="^88)
@printf("%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n",
        "z (m)", "Observed M[Ri_zz]", "Constitutive", "Coordinate", "Exact Total", "Error E_Δz")
println("-"^88)
for i in 1:n
    @printf("%6.1f | %17.6f | %12.6f | %12.6f | %12.6f | %12.6f\n",
            z_tower[i], M_Ri_zz[i], C_const[i], C_coord[i], Ri_exact[i], E_error[i])
end
println("="^88)