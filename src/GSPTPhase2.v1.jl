module GSPTPhase2

using LinearAlgebra, Statistics

export ProfileData, CoordinateGeometry, ConstitutiveGeometry, ObservationDiagnostic,
    DiagnosticResult, DomainMetrics, compute_gspt, compare_tracks

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
    sigma_obs::Float64         # Additive Ri-space observation noise std dev
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
    closure_residual::Vector{Float64}# Δ_closure = Ri_zz_obs - Ri_zz_gspt
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
    bias_R_scm::Vector{Float64} # Vertical profile of R_coord bias
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
    D1[n, (n-2):n] .= stencil_weights(z[(n-2):n], z[n], 2) # Boundary order
    return D1, D2
end

function solve_morozov(y_obs::Vector{Float64}, R_tilde::Matrix{Float64}, σ::Float64)
    n = length(y_obs)
    σ == 0.0 && return 0.0
    target = n * (σ^2)
    λ_grid = 10.0 .^ range(-8, 1, length=1000)
    best_λ = λ_grid[1]
    min_diff = Inf
    for λ in λ_grid
        y_s = (I(n) + λ .* R_tilde) \ y_obs
        res = sum((y_s .- y_obs) .^ 2)
        diff = abs(res - target)
        if diff < min_diff
            min_diff = diff
            best_λ = λ
        end
    end
    return best_λ
end

# --- 3. Independent Triple-Point Diagnostic ---
function compute_c_tp(z::Vector{Float64}, tke::Vector{Float64}, Km::Vector{Float64},
    D1::Matrix{Float64}; tke_fraction::Float64=0.05)
    H_SBL = maximum(z)

    # 1. TKE floor height z_e
    e_min = minimum(tke) + tke_fraction * (maximum(tke) - minimum(tke))
    idx_e = findfirst(e -> e <= e_min, tke)
    z_e = isnothing(idx_e) ? z[end] : z[idx_e]

    # 2. TKE gradient maximum height z_ez
    de_dz = abs.(D1 * tke)
    z_ez = z[argmax(de_dz)]

    # 3. Independent diffusivity extinction height z_K
    Km_min = minimum(Km) + tke_fraction * (maximum(Km) - minimum(Km))
    idx_K = findfirst(K -> K <= Km_min, Km)
    z_K = isnothing(idx_K) ? z[end] : z[idx_K]

    pts = [z_e, z_ez, z_K]
    return (maximum(pts) - minimum(pts)) / H_SBL
end

# --- 4. Core Computation Engine ---
function compute_gspt(data::ProfileData; β_m=1.0, β_h=3.0, eps_C=1e-5,
    is_observation=false, tke_fraction=0.05)
    z = data.z
    n = length(z)
    g, κ = 9.81, 0.40

    D1_dim, D2_dim = build_operators(z)

    # Obukhov Length L(z) using Virtual Potential Temperature
    ustar = ((data.uw .^ 2) .+ (data.vw .^ 2)) .^ 0.25
    L = zeros(n)
    for i in 1:n
        flux = data.wthv[i]
        L[i] = abs(flux) > 1e-6 ? -(ustar[i]^3 * data.theta_v[i]) / (κ * g * flux) : 1e5
    end

    # Coordinate Geometry Object
    zeta = z ./ L
    zeta_z = D1_dim * zeta
    zeta_zz = D2_dim * zeta
    coord_geom = CoordinateGeometry(zeta, zeta_z, zeta_zz)

    # Analytic Stability Parameterization (Ri_gspt)
    Ri_func(ζ) = ζ * (1.0 + β_h * ζ) / ((1.0 + β_m * ζ)^2)
    Ri_z_func(ζ) = (1.0 + (2.0 * β_h - β_m) * ζ) / ((1.0 + β_m * ζ)^3)
    Ri_zz_func(ζ) = (2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ) / ((1.0 + β_m * ζ)^4)

    Ri_gspt_val = Ri_func.(zeta)
    Ri_z_gspt_val = Ri_z_func.(zeta)
    Ri_zz_gspt_val = Ri_zz_func.(zeta)

    C_const = Ri_zz_gspt_val .* (zeta_z .^ 2)
    C_coord = Ri_z_gspt_val .* zeta_zz
    Ri_zz_chain_rule = C_const .+ C_coord

    R_coord = [abs(C_const[i]) > eps_C ? C_coord[i] / C_const[i] : NaN for i in 1:n]
    const_geom = ConstitutiveGeometry(Ri_gspt_val, Ri_z_gspt_val, Ri_zz_chain_rule, C_const, C_coord, R_coord)

    # Observed Diagnostic Object
    du_dz = D1_dim * data.u
    dv_dz = D1_dim * data.v
    dth_dz = D1_dim * data.theta_v
    S2 = max.((du_dz .^ 2) .+ (dv_dz .^ 2), 1e-6)
    Ri_raw = (g ./ data.theta_v) .* dth_dz ./ S2

    Ri_eval = zeros(n)
    if is_observation && data.sigma_obs > 0.0
        H_SBL = maximum(z)
        _, D2_tilde = build_operators(z ./ H_SBL)
        λ_opt = solve_morozov(Ri_raw, D2_tilde' * D2_tilde, data.sigma_obs)
        Ri_eval .= (I(n) + λ_opt .* (D2_tilde' * D2_tilde)) \ Ri_raw
    else
        Ri_eval .= Ri_raw
    end

    Ri_zz_obs = D2_dim * Ri_eval
    closure_residual = Ri_zz_obs .- Ri_zz_chain_rule
    tau_truncation = (D2_dim * Ri_gspt_val) .- Ri_zz_chain_rule

    obs_diag = ObservationDiagnostic(Ri_eval, Ri_zz_obs, tau_truncation, closure_residual)
    C_TP = compute_c_tp(z, data.tke, data.Km, D1_dim; tke_fraction=tke_fraction)

    return DiagnosticResult(z, L, coord_geom, const_geom, obs_diag, C_TP)
end

# --- 5. Multi-Track Comparison & Bias Evaluation ---
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

    # Common Valid Joint Masks
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