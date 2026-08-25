using LinearAlgebra, Printf, Random, Statistics

# --- 1. Non-Uniform Stencil Weights ---
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p);
    b[m+1] = 1.0
    return A \ b
end

function build_operators(z::Vector{Float64})
    n = length(z)
    D1, D2 = zeros(n, n), zeros(n, n)
    for i in 2:(n-1)
        idx = [i-1, i, i+1]
        D1[i, idx] .= stencil_weights(z[idx], z[i], 1)
        D2[i, idx] .= stencil_weights(z[idx], z[i], 2)
    end
    D1[1, 1:3] .= stencil_weights(z[1:3], z[1], 1)
    D2[1, 1:3] .= stencil_weights(z[1:3], z[1], 2)
    D1[n, (n-2):n] .= stencil_weights(z[(n-2):n], z[n], 1)
    D2[n, (n-2):n] .= stencil_weights(z[(n-2):n], z[n], 2)
    return D1, D2
end

# --- 2. Fully Analytical LLJ Kinematics & Stability Functions ---
β_m, β_h = 1.0, 3.0
L0, H_LLJ = 20.0, 60.0

L_func(z) = L0 * (1.0 - 0.45 * sin(π * z / H_LLJ))
L_p(z) = -L0 * 0.45 * (π / H_LLJ) * cos(π * z / H_LLJ)
L_pp(z) = L0 * 0.45 * (π / H_LLJ)^2 * sin(π * z / H_LLJ)

ζ_func(z) = z / L_func(z)
ζ_z_func(z) = (L_func(z) - z * L_p(z)) / (L_func(z)^2)
ζ_zz_func(z) = -2.0 * (L_p(z) / L_func(z)) * ζ_z_func(z) - (z * L_pp(z)) / (L_func(z)^2)

Ri_func(ζ) = ζ * (1.0 + β_h * ζ) / (1.0 + β_m * ζ)^2
Ri_z_func(ζ) = (1.0 + (2.0 * β_h - β_m) * ζ) / (1.0 + β_m * ζ)^3
Ri_zz_func(ζ) = (2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ) / (1.0 + β_m * ζ)^4

# --- 3. Tower Geometry & Nondimensional Setup ---
z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
n = length(z_tower)
H_SBL = 55.0
z_tilde = z_tower ./ H_SBL

D1_dim, D2_dim = build_operators(z_tower)
_, D2_tilde = build_operators(z_tilde)

R_tilde = D2_tilde' * D2_tilde

# Analytical Reference Profiles
ζ_true = ζ_func.(z_tower)
ζ_z_true = ζ_z_func.(z_tower)
ζ_zz_true = ζ_zz_func.(z_tower)

Ri_z_true = Ri_z_func.(ζ_true)
Ri_zz_true = Ri_zz_func.(ζ_true)

C_const_ana = Ri_zz_true .* (ζ_z_true .^ 2)
C_coord_ana = Ri_z_true .* ζ_zz_true
Ri_zz_ana = C_const_ana .+ C_coord_ana

# Diagnostic Ratio R_coord
R_coord = C_coord_ana ./ C_const_ana

# --- 4. Morozov Regularization (Dimensionless λ_tilde) ---
σ_obs = 0.008
Random.seed!(42)
Ri_true_vals = Ri_func.(ζ_true)
Ri_obs = Ri_true_vals .+ σ_obs .* randn(n)

function solve_morozov(y_obs, R_t, σ, n_pts)
    target = n_pts * (σ^2)
    λ_grid = 10.0 .^ range(-6, 2, length=1000)
    best_λ = λ_grid[1]
    min_diff = Inf
    for λ in λ_grid
        y_s = (I(n_pts) + λ .* R_t) \ y_obs
        res = sum((y_s .- y_obs) .^ 2)
        diff = abs(res - target)
        if diff < min_diff
            min_diff = diff
            best_λ = λ
        end
    end
    return best_λ
end

λ_tilde_opt = solve_morozov(Ri_obs, R_tilde, σ_obs, n)
Ri_smooth = (I(n) + λ_tilde_opt .* R_tilde) \ Ri_obs

# Observed Physical Curvature (Dimensional)
M_Ri_zz = D2_dim * Ri_smooth

# --- 5. Monte Carlo Resolution Error Study ---
N_mc = 2000
mc_errors = zeros(n, N_mc)

for k in 1:N_mc
    y_noisy = Ri_true_vals .+ σ_obs .* randn(n)
    λ_k = solve_morozov(y_noisy, R_tilde, σ_obs, n)
    y_sk = (I(n) + λ_k .* R_tilde) \ y_noisy
    mc_errors[:, k] .= (D2_dim * y_sk) .- Ri_zz_ana
end

rmse_z = sqrt.(mean(mc_errors .^ 2, dims=2))[:]

# --- 6. Print Verification Results ---
println("="^105)
@printf("%-6s | %-10s | %-12s | %-12s | %-12s | %-10s | %-10s | %-10s\n",
    "z (m)", "ζ_z", "C_const", "C_coord", "Ri_zz (Ana)", "R_coord", "M[Ri_zz]", "RMSE (MC)")
println("-"^105)
for i in 1:n
    @printf("%6.1f | %10.5f | %12.6f | %12.6f | %12.6f | %10.4f | %10.6f | %10.6f\n",
        z_tower[i], ζ_z_true[i], C_const_ana[i], C_coord_ana[i], Ri_zz_ana[i], R_coord[i], M_Ri_zz[i], rmse_z[i])
end
println("="^105)
@printf("Optimized Dimensionless Tikhonov Parameter (Morozov): λ_tilde = %.6f\n", λ_tilde_opt)