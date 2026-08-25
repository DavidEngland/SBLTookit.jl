module GSPTPhase2

using LinearAlgebra, Statistics

export ProfileData, CoordinateGeometry, ConstitutiveGeometry, ObservationDiagnostic,
    DiagnosticResult, DomainMetrics, compute_gspt, compare_tracks, check_tangential_cone

# --- 1. Data Structures ---
struct ProfileData
    z::Vector{Float64}         # Height array (m)
    u::Vector{Float64}         # Zonal velocity (m/s)
    v::Vector{Float64}         # Meridional velocity (m/s)
    theta_v::Vector{Float64}   # Virtual potential temperature (K)
    wthv::Vector{Float64}      # Kinematic virtual heat flux (K m/s)
    uw::Vector{Float64}        # Kinematic u-momentum flux (m^2/s^2)
    vw::Vector{Float64}        # Kinematic v-momentum flux (m^2/s^2)
    tke::Vector{Float64}       # Turbulent Kinetic Energy (m^2/s^2)
    Km::Vector{Float64}        # Turbulent diffusivity for momentum (m^2/s)
    sigma_u::Float64           # Sensor std dev for u (m/s)
    sigma_v::Float64           # Sensor std dev for v (m/s)
    sigma_theta::Float64       # Sensor std dev for theta_v (K)
    sigma_ri_diag::Float64     # Optional diagnostic-space Ri noise std dev
end

struct CoordinateGeometry
    zeta::Vector{Float64}
    zeta_z::Vector{Float64}
    zeta_zz::Vector{Float64}
end

struct ConstitutiveGeometry
    Ri_gspt::Vector{Float64}
    Ri_z_gspt::Vector{Float64}
    Ri_zz_gspt::Vector{Float64}
    C_const::Vector{Float64}
    C_coord::Vector{Float64}
    R_coord::Vector{Float64}    # Masked diagnostic (NaN where |C_const| <= eps_C)
end

struct ObservationDiagnostic
    Ri_obs::Vector{Float64}
    Ri_zz_obs::Vector{Float64}
    tau_truncation::Vector{Float64}  # Discrete Taylor truncation error proxy
    tau_reg_sens::Vector{Float64}    # Field-space vs Diagnostic-space MDP discrepancy
    closure_residual::Vector{Float64}# Δ_closure = Ri_zz_obs - Ri_zz_gspt
    ill_conditioned_mask::Vector{Bool} # True where S^2 <= S^2_min
    grid_cond_number::Float64       # κ(R_tilde) = λ_max / λ_min of discretization operator
end

struct DiagnosticResult
    z::Vector{Float64}
    L::Vector{Float64}
    coord_geom::CoordinateGeometry
    const_geom::ConstitutiveGeometry
    obs_diag::ObservationDiagnostic
    C_TP::Float64               # Triple-point spatial spread ratio
end

struct DomainMetrics
    masking_fraction_obs::Float64
    masking_fraction_scm::Float64
    masking_fraction_les::Float64
    bias_R_scm::Vector{Float64}
    bias_R_les::Vector{Float64}
    rmse_scm::Float64
    rmse_les::Float64
    mae_scm::Float64
    mae_les::Float64
    median_abs_diff_scm::Float64
    median_abs_diff_les::Float64
end

# --- 2. Stencil & Operator Utilities ---
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

function solve_morozov(y_obs::Vector{Float64}, R_tilde::Matrix{Float64}, σ::Float64)
    n = length(y_obs)
    σ == 0.0 && return copy(y_obs)
    target = n * (σ^2)

    # Stabilize matrix conditioning for inversion
    R_stable = R_tilde + 1e-8 * I(n)

    λ_grid = 10.0 .^ range(-8, 2, length=2000)
    best_λ = λ_grid[1]
    min_diff = Inf

    for λ in λ_grid
        y_s = (I(n) + λ .* R_stable) \ y_obs
        res = sum((y_s .- y_obs) .^ 2)
        diff = abs(res - target)
        if diff < min_diff
            min_diff = diff
            best_λ = λ
        end
    end
    return (I(n) + best_λ .* R_stable) \ y_obs
end

# --- 3. Tangential Cone Condition Evaluator ---
function check_tangential_cone(x1::Vector{Float64}, x2::Vector{Float64}, g::Float64, theta0::Float64, D1::Matrix{Float64})
    # Forward map F maps profile state x = [u; v; theta_v] to Ri profile
    n = length(x1) ÷ 3
    eval_F(x) = begin
        u, v, th = x[1:n], x[(n+1):(2*n)], x[(2*n+1):(3*n)]
        du = D1 * u;
        dv = D1 * v;
        dth = D1 * th
        S2 = max.((du .^ 2) .+ (dv .^ 2), 1e-6)
        return (g ./ theta0) .* dth ./ S2
    end

    F1 = eval_F(x1)
    F2 = eval_F(x2)

    # Numerical Jacobian calculation for F'(x1)
    eps_j = 1e-6
    J = zeros(n, 3*n)
    for j in 1:(3*n)
        x_plus = copy(x1);
        x_plus[j] += eps_j
        J[:, j] = (eval_F(x_plus) .- F1) ./ eps_j
    end

    dx = x2 .- x1
    linear_approx = F1 .+ J * dx

    denom = norm(F2 .- F1)
    denom == 0.0 && return 0.0
    gamma = norm(F2 .- linear_approx) / denom
    return gamma
end

# --- 4. Independent Triple-Point Diagnostic ---
function compute_c_tp(z::Vector{Float64}, tke::Vector{Float64}, Km::Vector{Float64},
    D1::Matrix{Float64}; tke_fraction::Float64=0.05)
    H_SBL = maximum(z)

    e_min = minimum(tke) + tke_fraction * (maximum(tke) - minimum(tke))
    idx_e = findfirst(e -> e <= e_min, tke)
    z_e = isnothing(idx_e) ? z[end] : z[idx_e]

    de_dz = abs.(D1 * tke)
    z_ez = z[argmax(de_dz)]

    Km_min = minimum(Km) + tke_fraction * (maximum(Km) - minimum(Km))
    idx_K = findfirst(K -> K <= Km_min, Km)
    z_K = isnothing(idx_K) ? z[end] : z[idx_K]

    pts = [z_e, z_ez, z_K]
    return (maximum(pts) - minimum(pts)) / H_SBL
end

# --- 5. Core Computation Engine ---
function compute_gspt(data::ProfileData; β_m=1.0, β_h=3.0, eps_C=1e-5, S2_min=1e-4,
    is_observation=false, tke_fraction=0.05, mask_ill_conditioned::Bool=true)
    z = data.z
    n = length(z)
    g, κ = 9.81, 0.40
    H_SBL = maximum(z)

    D1_dim, D2_dim = build_operators(z)
    _, D2_tilde = build_operators(z ./ H_SBL)
    R_tilde = D2_tilde' * D2_tilde

    # Calculate grid operator condition number κ(R_tilde)
    eigs = eigvals(R_tilde)
    valid_eigs = filter(e -> e > 1e-12, eigs)
    grid_cond_number = length(valid_eigs) > 0 ? maximum(valid_eigs) / minimum(valid_eigs) : Inf

    # Track A: Primitive Field Regularization (Primary MDP Track)
    u_eval = is_observation ? solve_morozov(data.u, R_tilde, data.sigma_u) : copy(data.u)
    v_eval = is_observation ? solve_morozov(data.v, R_tilde, data.sigma_v) : copy(data.v)
    th_eval = is_observation ? solve_morozov(data.theta_v, R_tilde, data.sigma_theta) : copy(data.theta_v)

    du_dz = D1_dim * u_eval
    dv_dz = D1_dim * v_eval
    dth_dz = D1_dim * th_eval

    S2 = (du_dz .^ 2) .+ (dv_dz .^ 2)
    ill_conditioned_mask = S2 .<= S2_min
    S2_bounded = max.(S2, S2_min)

    Ri_obs_field = (g ./ th_eval) .* dth_dz ./ S2_bounded
    Ri_zz_obs_field = D2_dim * Ri_obs_field

    # Track B: Diagnostic-Space Regularization (Sensitivity Benchmark)
    du_raw = D1_dim * data.u;
    dv_raw = D1_dim * data.v;
    dth_raw = D1_dim * data.theta_v
    S2_raw = max.((du_raw .^ 2) .+ (dv_raw .^ 2), S2_min)
    Ri_raw = (g ./ data.theta_v) .* dth_raw ./ S2_raw
    Ri_obs_diag = is_observation && data.sigma_ri_diag > 0.0 ? solve_morozov(Ri_raw, R_tilde, data.sigma_ri_diag) : copy(Ri_raw)
    Ri_zz_obs_diag = D2_dim * Ri_obs_diag

    tau_reg_sens = Ri_zz_obs_field .- Ri_zz_obs_diag

    # Obukhov Length L(z)
    ustar = ((data.uw .^ 2) .+ (data.vw .^ 2)) .^ 0.25
    L = zeros(n)
    for i in 1:n
        flux = data.wthv[i]
        L[i] = abs(flux) > 1e-6 ? -(ustar[i]^3 * th_eval[i]) / (κ * g * flux) : 1e5
    end

    # Coordinate Geometry
    zeta = z ./ L
    zeta_z = D1_dim * zeta
    zeta_zz = D2_dim * zeta
    coord_geom = CoordinateGeometry(zeta, zeta_z, zeta_zz)

    # Analytic Parameterization
    Ri_func(ζ) = ζ * (1.0 + β_h * ζ) / ((1.0 + β_m * ζ)^2)
    Ri_z_func(ζ) = (1.0 + (2.0 * β_h - β_m) * ζ) / ((1.0 + β_m * ζ)^3)
    Ri_zz_func(ζ) = (2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ) / ((1.0 + β_m * ζ)^4)

    Ri_gspt_val = Ri_func.(zeta)
    Ri_z_gspt_val = Ri_z_func.(zeta)
    Ri_zz_gspt_val = Ri_zz_func.(zeta)

    C_const = Ri_zz_gspt_val .* (zeta_z .^ 2)
    C_coord = Ri_z_gspt_val .* zeta_zz
    Ri_zz_chain_rule = C_const .+ C_coord

    R_coord = [
        abs(C_const[i]) > eps_C && (!mask_ill_conditioned || !ill_conditioned_mask[i]) ? C_coord[i] / C_const[i] : NaN
        for i in 1:n
    ]
    const_geom = ConstitutiveGeometry(Ri_gspt_val, Ri_z_gspt_val, Ri_zz_chain_rule, C_const, C_coord, R_coord)

    closure_residual = Ri_zz_obs_field .- Ri_zz_chain_rule
    tau_truncation = (D2_dim * Ri_gspt_val) .- Ri_zz_chain_rule

    obs_diag = ObservationDiagnostic(Ri_obs_field, Ri_zz_obs_field, tau_truncation,
        tau_reg_sens, closure_residual, ill_conditioned_mask, grid_cond_number)

    C_TP = compute_c_tp(z, data.tke, data.Km, D1_dim; tke_fraction=tke_fraction)

    return DiagnosticResult(z, L, coord_geom, const_geom, obs_diag, C_TP)
end

# --- 6. Multi-Track Comparison ---
function compare_tracks(obs::DiagnosticResult, scm::DiagnosticResult, les::DiagnosticResult;
    eps_mask=0.15)
    n = length(obs.z)

    calc_mask_frac(R_vec) = begin
        v = filter(!isnan, R_vec)
        length(v) > 0 ? count(r -> abs(r + 1.0) < eps_mask, v) / length(v) : 0.0
    end

    frac_obs = calc_mask_frac(obs.const_geom.R_coord)
    frac_scm = calc_mask_frac(scm.const_geom.R_coord)
    frac_les = calc_mask_frac(les.const_geom.R_coord)

    mask_scm = .!isnan.(obs.const_geom.R_coord) .& .!isnan.(scm.const_geom.R_coord)
    mask_les = .!isnan.(obs.const_geom.R_coord) .& .!isnan.(les.const_geom.R_coord)

    bias_R_scm = fill(NaN, n)
    bias_R_les = fill(NaN, n)

    bias_R_scm[mask_scm] .= scm.const_geom.R_coord[mask_scm] .- obs.const_geom.R_coord[mask_scm]
    bias_R_les[mask_les] .= les.const_geom.R_coord[mask_les] .- obs.const_geom.R_coord[mask_les]

    diff_scm = filter(!isnan, bias_R_scm)
    diff_les = filter(!isnan, bias_R_les)

    rmse_scm = sqrt(mean(diff_scm .^ 2))
    rmse_les = sqrt(mean(diff_les .^ 2))
    mae_scm = mean(abs.(diff_scm))
    mae_les = mean(abs.(diff_les))
    med_scm = median(abs.(diff_scm))
    med_les = median(abs.(diff_les))

    return DomainMetrics(frac_obs, frac_scm, frac_les, bias_R_scm, bias_R_les,
        rmse_scm, rmse_les, mae_scm, mae_les, med_scm, med_les)
end

end # module