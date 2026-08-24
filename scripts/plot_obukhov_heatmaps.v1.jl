#!/usr/bin/env julia
# SBLToolkit.jl/scripts/plot_obukhov_heatmaps.jl

using Pkg

# -------------------------------------------------------------------
# 1. Environment & Path Resolution (Directory Independent)
# -------------------------------------------------------------------
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

using SBLToolkit
using CSV
using DataFrames
using Dates
using Plots
using Plots.PlotMeasures
using Statistics

const CAMPAIGN_ROOT = normpath(joinpath(PROJECT_ROOT, "..", "SpectralBL-Analytics"))
const CAMPAIGN_SOURCES = Dict(
    "CASES-99" => joinpath(CAMPAIGN_ROOT, "data", "drafts", "trajectories", "trajectory_cases_99.csv"),
    "FLOSS" => joinpath(CAMPAIGN_ROOT, "data", "drafts", "trajectories", "trajectory_floss.csv"),
    "BLLAST" => joinpath(CAMPAIGN_ROOT, "data", "drafts", "trajectories", "trajectory_bllast.csv"),
    "SHEBA" => joinpath(CAMPAIGN_ROOT, "data", "drafts", "trajectories", "trajectory_sheba.csv"),
)

const SBL_OUTPUT_DIR = joinpath(PROJECT_ROOT, "reports", "generated", "sbltoolkit_heatmaps")
mkpath(SBL_OUTPUT_DIR)

slugify(s::String) = lowercase(replace(s, r"[^A-Za-z0-9]+" => "_"))

function to_f64(x)
    (x === missing || x === nothing) && return NaN
    x isa Number && return Float64(x)
    s = strip(String(x))
    (isempty(s) || lowercase(s) == "nan") && return NaN
    v = tryparse(Float64, s)
    return v === nothing ? NaN : v
end

# -------------------------------------------------------------------
# 2. Robust Non-Uniform Stencils (NaN-Mask Preserving)
# -------------------------------------------------------------------
function non_uniform_gradient_1d(f::Vector{Float64}, z::Vector{Float64})
    n = length(z)
    df = fill(NaN, n)
    valid_idx = findall(isfinite, f)
    length(valid_idx) < 3 && return df

    zv = z[valid_idx]
    fv = f[valid_idx]
    nv = length(zv)
    dfv = zeros(Float64, nv)

    # 3-point non-uniform finite difference stencil
    dfv[1] = (fv[2] - fv[1]) / (zv[2] - zv[1])
    dfv[end] = (fv[end] - fv[end-1]) / (zv[end] - zv[end-1])

    for i in 2:(nv-1)
        h1 = zv[i] - zv[i-1]
        h2 = zv[i+1] - zv[i]
        dfv[i] = (fv[i+1]*h1^2 - fv[i-1]*h2^2 + fv[i]*(h2^2 - h1^2)) / (h1 * h2 * (h1 + h2))
    end

    df[valid_idx] .= dfv
    return df
end

function non_uniform_hessian_1d(f::Vector{Float64}, z::Vector{Float64})
    df1 = non_uniform_gradient_1d(f, z)
    return non_uniform_gradient_1d(df1, z)
end

# -------------------------------------------------------------------
# 3. Data Ingestion & Physical Matrix Assembly
# -------------------------------------------------------------------
function extract_obukhov_matrix(df::DataFrame; g=9.81, theta_0=273.15, k_vk=0.40)
    time_col = hasproperty(df, :time) ? :time : (hasproperty(df, :sample_index) ? :sample_index : nothing)
    z_col = hasproperty(df, :z) ? :z : (hasproperty(df, :height) ? :height : nothing)

    @assert time_col !== nothing "DataFrame requires a valid time or sample_index column."
    @assert z_col !== nothing "No vertical coordinate found in dataset; cannot fabricate tower heights."

    timestamps = sort(unique(df[!, time_col]))
    z_levels = sort(unique(filter(isfinite, [to_f64(v) for v in df[!, z_col]])))

    @assert !isempty(z_levels) "Vertical coordinate column contains no valid numeric values."

    nz, nt = length(z_levels), length(timestamps)
    L_mat = fill(NaN, nz, nt)

    time_map = Dict(t => j for (j, t) in enumerate(timestamps))
    z_map = Dict(z => i for (i, z) in enumerate(z_levels))

    l_sym = hasproperty(df, :L_obukhov) ? :L_obukhov : (hasproperty(df, :L) ? :L : nothing)

    if l_sym !== nothing
        for row in eachrow(df)
            t_val = row[time_col]
            z_val = to_f64(row[z_col])
            (haskey(time_map, t_val) && haskey(z_map, z_val)) || continue
            L_mat[z_map[z_val], time_map[t_val]] = to_f64(row[l_sym])
        end
    elseif all(c -> hasproperty(df, c), [:u_star, :sensible_heat_flux])
        # Kinematic Heat Flux Conversion: H [W/m²] / (ρ * c_p) where ρ*c_p ≈ 1216.0 J/(m³ K)
        rho_cp = 1216.0
        for row in eachrow(df)
            t_val = row[time_col]
            z_val = to_f64(row[z_col])
            (haskey(time_map, t_val) && haskey(z_map, z_val)) || continue

            ustar = to_f64(row[:u_star])
            shf = to_f64(row[:sensible_heat_flux])

            if isfinite(ustar) && isfinite(shf) && abs(shf) > 1e-4
                w_theta = shf / rho_cp
                L_mat[z_map[z_val], time_map[t_val]] = -(ustar^3 * theta_0) / (k_vk * g * w_theta)
            end
        end
    end

    # Preserves missingness as NaN; converts valid finite L to 1/L
    inv_L_mat = fill(NaN, nz, nt)
    for j in 1:nt, i in 1:nz
        val = L_mat[i, j]
        if isfinite(val)
            inv_L_mat[i, j] = abs(val) > 1e-6 ? 1.0 / val : 0.0 # Neutral asymptote
        end
    end

    # Calculate ζ = z / L = z * (1/L), preserving NaNs automatically
    zeta_mat = fill(NaN, nz, nt)
    for j in 1:nt, i in 1:nz
        if isfinite(inv_L_mat[i, j])
            zeta_mat[i, j] = z_levels[i] * inv_L_mat[i, j]
        end
    end

    # Calculate coordinate Jacobian ζ_z = ∂ζ/∂z and curvature ζ_zz = ∂²ζ/∂z²
    zeta_z_mat = fill(NaN, nz, nt)
    zeta_zz_mat = fill(NaN, nz, nt)
    for j in 1:nt
        zeta_z_mat[:, j] = non_uniform_gradient_1d(zeta_mat[:, j], z_levels)
        zeta_zz_mat[:, j] = non_uniform_hessian_1d(zeta_mat[:, j], z_levels)
    end

    return (
        timestamps=timestamps,
        z_levels=z_levels,
        inv_L=inv_L_mat,
        zeta=zeta_mat,
        zeta_z=zeta_z_mat,
        zeta_zz=zeta_zz_mat
    )
end

# -------------------------------------------------------------------
# 4. Geometric Manifold Diagnostic Plotting
# -------------------------------------------------------------------
function plot_sbltoolkit_obukhov_panel(data::NamedTuple, campaign_name::String, out_path::String)
    t_axis = data.timestamps
    z_axis = data.z_levels

    opts = (
        ylabel="Height z [m]",
        framestyle=:box,
        guidefontsize=9,
        tickfontsize=8,
        titlefontsize=10,
        margin=3mm
    )

    # Panel 1: Signed Inverse Obukhov Scale 1/L(z,t)
    p1 = heatmap(t_axis, z_axis, clamp.(data.inv_L, -0.1, 0.1);
        clims=(-0.05, 0.05), color=:balance,
        title="$(campaign_name) — Inverse Scale: 1/L(z,t)",
        colorbar_title=" 1/L [m⁻¹]", opts...)

    # Panel 2: Dimensionless Stability Coordinate ζ(z,t) = z/L
    zeta_symlog = @. sign(data.zeta) * log10(1.0 + abs(data.zeta))
    p2 = heatmap(t_axis, z_axis, zeta_symlog;
        clims=(-1.5, 1.5), color=:pu_or,
        title="$(campaign_name) — Stability Profile: sgn(ζ) log₁₀(1 + |ζ|)",
        colorbar_title=" sgn(ζ) log₁₀(1+|ζ|)", opts...)

    # Panel 3: Coordinate Jacobian ζ_z = ∂(z/L)/∂z with Fold Contour (ζ_z = 0)
    p3 = heatmap(t_axis, z_axis, clamp.(data.zeta_z, -0.2, 0.2);
        clims=(-0.1, 0.1), color=:vik,
        title="$(campaign_name) — Jacobian: ζ_z = ∂(z/L)/∂z",
        colorbar_title=" ζ_z [m⁻¹]", opts...)

    # Overlay fold line condition ζ_z = 0
    if count(isfinite, data.zeta_z) > 0
        contour!(p3, t_axis, z_axis, data.zeta_z, levels=[0.0], color=:white, lw=1.5, ls=:solid)
    end

    # Panel 4: Coordinate Curvature ζ_zz = ∂²(z/L)/∂z²
    p4 = heatmap(t_axis, z_axis, clamp.(data.zeta_zz, -0.05, 0.05);
        clims=(-0.02, 0.02), color=:curl,
        title="$(campaign_name) — Curvature: ζ_zz = ∂²(z/L)/∂z²",
        colorbar_title=" ζ_zz [m⁻²]", opts...)

    fig = plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 750), dpi=300)
    savefig(fig, out_path)
    return fig
end

# -------------------------------------------------------------------
# 5. Execution Loop
# -------------------------------------------------------------------
function generate_all_sbltoolkit_heatmaps()
    for (name, path) in CAMPAIGN_SOURCES
        if !isfile(path)
            @warn "Campaign file not found: $path. Skipping."
            continue
        end

        df = CSV.read(path, DataFrame)
        obukhov_data = extract_obukhov_matrix(df)

        fig_out = joinpath(SBL_OUTPUT_DIR, "$(slugify(name))_obukhov_heatmaps.png")
        plot_sbltoolkit_obukhov_panel(obukhov_data, name, fig_out)
        @info "Rendered geometric heatmaps for $(name) -> $fig_out"
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_all_sbltoolkit_heatmaps()
end