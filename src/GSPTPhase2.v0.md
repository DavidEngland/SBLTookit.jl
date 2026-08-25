### GABLS3 Phase 2 Architecture: Comparison Layer

```
                        GABLS3 Ingestion Layer
                    (time × height state fields)
                               │
                       Matched Vertical Grid
                               │
               ┌───────────────┴───────────────┐
               ▼                               ▼
       Observation Track                  Model Tracks
   (Regularized: M_Δz via λ)        (Unregularized: τ_Δz)
               │                               │
               └───────────────┬───────────────┘
                               ▼
               ┌───────────────────────────────┐
               │    GSPT Core Engine           │
               │  - L(z), ζ, ζ_z, ζ_zz          │
               │  - Ri, Ri_ζ, Ri_ζζ             │
               │  - C_const, C_coord            │
               │  - Masked R_coord              │
               │  - Triple-Point C_TP           │
               └───────────────┬───────────────┘
                               ▼
               ┌───────────────────────────────┐
               │ Multi-Track Comparison & Bias │
               │  - Footprints: R_coord(z,t)   │
               │  - Bias: B_R(z,t)             │
               │  - Masking Fraction (|R+1|<ε) │
               └───────────────────────────────┘

```

---

### Phase 2 Modular Julia Implementation

```julia
module GSPTPhase2

using LinearAlgebra, Statistics

export ProfileData, DiagnosticResult, compute_gspt, compare_tracks, DomainMetrics

# --- Data Containers ---
struct ProfileData
    z::Vector{Float64}         # Height array (m)
    u::Vector{Float64}         # Zonal velocity (m/s)
    v::Vector{Float64}         # Meridional velocity (m/s)
    theta::Vector{Float64}     # Potential temperature (K)
    wth::Vector{Float64}       # Kinematic heat flux (K m/s)
    uw::Vector{Float64}        # Kinematic u-momentum flux (m^2/s^2)
    vw::Vector{Float64}        # Kinematic v-momentum flux (m^2/s^2)
    tke::Vector{Float64}       # Turbulent Kinetic Energy (m^2/s^2)
    sigma_obs::Float64         # Instrument noise standard deviation (0.0 for models)
end

struct DiagnosticResult
    z::Vector{Float64}
    L::Vector{Float64}
    zeta::Vector{Float64}
    zeta_z::Vector{Float64}
    zeta_zz::Vector{Float64}
    Ri::Vector{Float64}
    C_const::Vector{Float64}
    C_coord::Vector{Float64}
    R_coord::Vector{Float64}    # Masked diagnostic (NaN where |C_const| <= eps_C)
    Ri_zz::Vector{Float64}
    C_TP::Float64               # Triple-point spatial spread ratio
end

struct DomainMetrics
    masking_fraction::Float64   # Fraction of points where |R_coord + 1| < eps_mask
    mean_bias_scm::Vector{Float64}
    mean_bias_les::Vector{Float64}
end

# --- Stencil & Operator Utilities ---
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p); b[m + 1] = 1.0
    return A \ b
end

function build_operators(z::Vector{Float64})
    n = length(z)
    D1, D2 = zeros(n, n), zeros(n, n)
    for i in 2:n-1
        idx = [i-1, i, i+1]
        D1[i, idx] .= stencil_weights(z[idx], z[i], 1)
        D2[i, idx] .= stencil_weights(z[idx], z[i], 2)
    end
    D1[1, 1:3] .= stencil_weights(z[1:3], z[1], 1)
    D2[1, 1:3] .= stencil_weights(z[1:3], z[1], 2)
    D1[n, n-2:n] .= stencil_weights(z[n-2:n], z[n], 1)
    D2[n, n-2:n] .= stencil_weights(z[n-2:n], z[n], 2)
    return D1, D2
end

function solve_morozov(y_obs::Vector{Float64}, R_tilde::Matrix{Float64}, σ::Float64)
    n = length(y_obs)
    σ == 0.0 && return 0.0 # No regularization for noise-free model data
    target = n * (σ^2)
    λ_grid = 10.0 .^ range(-8, 1, length=1000)
    best_λ = λ_grid[1]
    min_diff = Inf
    for λ in λ_grid
        y_s = (I(n) + λ .* R_tilde) \ y_obs
        res = sum((y_s .- y_obs).^2)
        diff = abs(res - target)
        if diff < min_diff
            min_diff = diff
            best_λ = λ
        end
    end
    return best_λ
end

# --- Triple-Point Convergence Diagnostic ---
function compute_c_tp(z::Vector{Float64}, tke::Vector{Float64}, D1::Matrix{Float64})
    H_SBL = maximum(z)
    n = length(z)

    # 1. TKE floor height z_e
    e_min = minimum(tke) + 0.05 * (maximum(tke) - minimum(tke))
    idx_e = findfirst(e -> e <= e_min, tke)
    z_e = isnothing(idx_e) ? z[end] : z[idx_e]

    # 2. TKE gradient maximum height z_ez
    de_dz = abs.(D1 * tke)
    idx_ez = argmax(de_dz)
    z_ez = z[idx_ez]

    # 3. Diffusivity extinction proxy z_K
    z_K = (z_e + z_ez) / 2.0

    pts = [z_e, z_ez, z_K]
    return (maximum(pts) - minimum(pts)) / H_SBL
end

# --- Core Agnostic GSPT Calculation Engine ---
function compute_gspt(data::ProfileData; β_m=1.0, β_h=3.0, eps_C=1e-5, is_observation=false)
    z = data.z
    n = length(z)
    g, κ = 9.81, 0.40

    D1_dim, D2_dim = build_operators(z)

    # Local Obukhov Length Profile
    ustar = ((data.uw .^ 2) .+ (data.vw .^ 2)) .^ 0.25
    L = zeros(n)
    for i in 1:n
        flux = data.wth[i]
        L[i] = abs(flux) > 1e-6 ? -(ustar[i]^3 * data.theta[i]) / (κ * g * flux) : 1e5
    end

    # Similarity Coordinate & Derivatives
    zeta = z ./ L
    zeta_z = D1_dim * zeta
    zeta_zz = D2_dim * zeta

    # Raw Richardson Number Profile
    du_dz = D1_dim * data.u
    dv_dz = D1_dim * data.v
    dth_dz = D1_dim * data.theta
    S2 = max.((du_dz.^2) .+ (dv_dz.^2), 1e-6)
    Ri_raw = (g ./ data.theta) .* dth_dz ./ S2

    # Dual-Track Processing: Noise Regularization vs Numerical Discretization
    Ri_eval = zeros(n)
    if is_observation && data.sigma_obs > 0.0
        H_SBL = maximum(z)
        _, D2_tilde = build_operators(z ./ H_SBL)
        R_tilde = D2_tilde' * D2_tilde
        λ_opt = solve_morozov(Ri_raw, R_tilde, data.sigma_obs)
        Ri_eval .= (I(n) + λ_opt .* R_tilde) \ Ri_raw
    else
        Ri_eval .= Ri_raw # Retain raw discretization bias τ_Δz for models
    end

    # Constitutive & Coordinate Curvature Components
    # Using analytic derivative formulations of Ri(ζ)
    Ri_z_func(ζ) = (1.0 + (2.0 * β_h - β_m) * ζ) / ((1.0 + β_m * ζ)^3)
    Ri_zz_func(ζ) = (2.0 * β_h - 4.0 * β_m + 2.0 * β_m * (β_m - 2.0 * β_h) * ζ) / ((1.0 + β_m * ζ)^4)

    Ri_z_val = Ri_z_func.(zeta)
    Ri_zz_val = Ri_zz_func.(zeta)

    C_const = Ri_zz_val .* (zeta_z .^ 2)
    C_coord = Ri_z_val .* zeta_zz
    Ri_zz_total = D2_dim * Ri_eval

    # Masked Diagnostic R_coord
    R_coord = zeros(n)
    for i in 1:n
        if abs(C_const[i]) > eps_C
            R_coord[i] = C_coord[i] / C_const[i]
        else
            R_coord[i] = NaN # Mask mathematical singularity
        end
    end

    C_TP = compute_c_tp(z, data.tke, D1_dim)

    return DiagnosticResult(z, L, zeta, zeta_z, zeta_zz, Ri_eval, C_const, C_coord, R_coord, Ri_zz_total, C_TP)
end

# --- Multi-Track Comparison Layer ---
function compare_tracks(obs::DiagnosticResult, scm::DiagnosticResult, les::DiagnosticResult; eps_mask=0.15)
    n = length(obs.z)

    # Direct Model Biases: B_R(z) = R_model - R_obs
    B_R_scm = scm.R_coord .- obs.R_coord
    B_R_les = les.R_coord .- obs.R_coord

    # Inflection Masking Test (|R_coord + 1| < eps_mask)
    valid_obs = filter(!isnan, obs.R_coord)
    count_masked = count(r -> abs(r + 1.0) < eps_mask, valid_obs)
    frac_masked = length(valid_obs) > 0 ? count_masked / length(valid_obs) : 0.0

    return DomainMetrics(frac_masked, B_R_scm, B_R_les)
end

end # module

```

---

### Implementation & Diagnostic Verification

| Target Metric | Formula | Mathematical Safeguard | Diagnostic Objective |
| --- | --- | --- | --- |
| **Masked Ratio ($\mathcal{R}_{\text{coord}}$)** | $\frac{C_{\text{coord}}}{C_{\text{const}}}$ | Evaluate as `NaN` if $\vert{}C_{\text{const}}\vert{} \le \epsilon_C$ | Eliminates division-by-zero singularities near simple inflection points. |
| **Direct Model Bias ($B_{\mathcal{R}}$)** | $\mathcal{R}_{\text{coord}}^{\text{model}} - \mathcal{R}_{\text{coord}}^{\text{obs}}$ | Computed on matched vertical grids | Quantifies misrepresentation of flux geometry across model runs. |
| **Masking Fraction** | $\frac{\text{count}(\vert{}\mathcal{R}_{\text{coord}} + 1\vert{} < \epsilon)}{\text{total valid points}}$ | Excludes `NaN` values from denominator | Measures domain extent controlled by curvature cancellation. |
| **Triple-Point Spread ($C_{\text{TP}}$)** | $\frac{\max\vert{}z_i - z_j\vert{}}{H_{\text{SBL}}}$ | Bound to domain boundary-layer scale | Distinguishes physical turbulent collapse ($C_{\text{TP}} \to 0$) from model numerical diffusion. |

---

### Key Scientific Output Tests

* **Inflection Masking Benchmark ($\mathcal{R}_{\text{coord}}^{\text{obs}} \approx -1.0 \land \mathcal{R}_{\text{coord}}^{\text{SCM}} \approx 0$):**
If observations yield $\mathcal{R}_{\text{coord}} \approx -1.0$ while the SCM reports $\mathcal{R}_{\text{coord}} \approx 0$, the SCM has preserved the constitutive stability function $\phi_h(\zeta)$ but underrepresented flux divergence ($L'(z)$). Retuning empirical coefficients $\beta_m, \beta_h$ in the model is unjustified; the parameterization instead requires non-local flux divergence acceleration terms.
* **Coordinate Fold Singularity Diagnostic ($\mathcal{R}_{\text{coord}} \gg 1.0$):**
Near the Low-Level Jet nose ($z \approx 80\text{--}140\text{ m}$), coordinate derivative $\zeta_z \to 0$ causes $C_{\text{const}} \to 0$ while $C_{\text{coord}}$ remains finite. The masked operator isolates this layer, preventing the model evaluation suite from interpreting profile flattening near the jet nose as turbulent decay.