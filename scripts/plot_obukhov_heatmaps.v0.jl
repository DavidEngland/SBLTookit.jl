#!/usr/bin/env julia
# SBLTookit.jl/scripts/plot_obukhov_heatmaps.jl

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

# Push local src/ into LOAD_PATH so Julia resolves the root package
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))

using SBLTookit
using CSV
using DataFrames
using Dates
using Plots
using Plots.PlotMeasures
using Statistics

# -------------------------------------------------------------------
# Campaign Data Source Paths & Output Directory Setup
# -------------------------------------------------------------------
const CAMPAIGN_SOURCES = Dict(
    "CASES-99" => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_cases_99.csv",
    "FLOSS"    => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_floss.csv",
    "BLLAST"   => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_bllast.csv",
    "SHEBA"    => "../SpectralBL-Analytics/data/drafts/trajectories/trajectory_sheba.csv",
)

const SBL_OUTPUT_DIR = joinpath(pwd(), "reports", "generated", "sbltookit_heatmaps")
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
# 1. Obukhov Matrix Assembly & Singular Transform
# -------------------------------------------------------------------
"""
    extract_obukhov_matrix(df::DataFrame)

Extracts or computes the 2D height-time Obukhov matrix L(z,t) and the
continuous curvature coordinate 1/L(z,t) across non-uniform tower levels.
"""
function extract_obukhov_matrix(df::DataFrame; g=9.81, theta_0=273.15, k_vk=0.40)
    time_col = hasproperty(df, :time) ? :time : (hasproperty(df, :sample_index) ? :sample_index : nothing)
    z_col = hasproperty(df, :z) ? :z : (hasproperty(df, :height) ? :height : nothing)

    @assert time_col !== nothing "DataFrame requires a valid time or sample_index column."

    timestamps = sort(unique(df[!, time_col]))

    # Identify vertical tower levels
    z_levels = if z_col !== nothing
        sort(unique(filter(isfinite, [to_f64(v) for v in df[!, z_col]])))
    else
        [1.5, 3.0, 6.0, 10.0, 20.0, 30.0, 50.0] # Fallback tower levels [m]
    end

    nz, nt = length(z_levels), length(timestamps)
    L_mat = fill(NaN, nz, nt)

    time_map = Dict(t => j for (j, t) in enumerate(timestamps))
    z_map = Dict(z => i for (i, z) in enumerate(z_levels))

    # Ingest direct L observations or flux profiles
    if hasproperty(df, :L_obukhov) || hasproperty(df, :L)
        l_sym = hasproperty(df, :L_obukhov) ? :L_obukhov : :L
        for row in eachrow(df)
            t_val = row[time_col]
            z_val = z_col !== nothing ? to_f64(row[z_col]) : z_levels[1]
            (haskey(time_map, t_val) && haskey(z_map, z_val)) || continue
            L_mat[z_map[z_val], time_map[t_val]] = to_f64(row[l_sym])
        end
    elseif all(c -> hasproperty(df, c), [:u_star, :sensible_heat_flux])
        for row in eachrow(df)
            t_val = row[time_col]
            z_val = z_col !== nothing ? to_f64(row[z_col]) : z_levels[1]
            (haskey(time_map, t_val) && haskey(z_map, z_val)) || continue

            ustar = to_f64(row[:u_star])
            shf = to_f64(row[:sensible_heat_flux])

            if isfinite(ustar) && isfinite(shf) && abs(shf) > 1e-4
                L_mat[z_map[z_val], time_map[t_val]] = -(ustar^3 * theta_0) / (k_vk * g * (shf / 1200.0))
            end
        end
    end

    # Compute inverse Obukhov length 1/L to eliminate division-by-zero singularities at neutral transitions
    inv_L_mat = @. 1.0 / L_mat
    inv_L_mat[isnan.(inv_L_mat)] .= 0.0 # Bounded representation for near-neutral L -> ∞

    # Log-transformed absolute Obukhov scale: sgn(L) * log10(|L|)
    log_L_mat = zeros(Float64, nz, nt)
    for j in 1:nt, i in 1:nz
        val = L_mat[i, j]
        if isfinite(val) && abs(val) > 0.1
            log_L_mat[i, j] = sign(val) * log10(abs(val))
        else
            log_L_mat[i, j] = NaN
        end
    end

    return (
        timestamps = timestamps,
        z_levels = z_levels,
        L = L_mat,
        inv_L = inv_L_mat,
        log_L = log_L_mat
    )
end

# -------------------------------------------------------------------
# 2. Comparative Campaign Obukhov Heat Map Panel Plotter
# -------------------------------------------------------------------
"""
    plot_sbltookit_obukhov_panel(data::NamedTuple, campaign_name::String, out_path::String)

Generates comparative 2x2 heatmaps displaying observational Obukhov scales,
inverse length 1/L(z,t), stability parameter ζ(z,t), and local profile divergence.
"""
function plot_sbltookit_obukhov_panel(data::NamedTuple, campaign_name::String, out_path::String)
    t_axis = 1:length(data.timestamps)
    z_axis = data.z_levels

    # Local stability parameter ζ(z,t) = z / L(z,t)
    nz, nt = length(z_axis), length(t_axis)
    zeta_mat = zeros(Float64, nz, nt)
    for j in 1:nt, i in 1:nz
        zeta_mat[i, j] = z_axis[i] * data.inv_L[i, j]
    end
    zeta_symlog = @. sign(zeta_mat) * log10(1.0 + abs(zeta_mat))

    opts = (
        ylabel = "Height z [m]",
        framestyle = :box,
        guidefontsize = 9,
        tickfontsize = 8,
        titlefontsize = 10,
        margin = 3mm
    )

    # Panel 1: Log-Scaled Obukhov Length Scale sgn(L) log₁₀|L|
    p1 = heatmap(t_axis, z_axis, data.log_L;
        clims = (-3.0, 3.0), color = :coolwarm,
        title = "$(campaign_name) — Obukhov Scale: sgn(L) log₁₀(|L|)",
        colorbar_title = " sgn(L) log₁₀|L| [m]", opts...)

    # Panel 2: Continuous Inverse Obukhov Curvature 1/L(z,t)
    p2 = heatmap(t_axis, z_axis, clamp.(data.inv_L, -0.1, 0.1);
        clims = (-0.05, 0.05), color = :balance,
        title = "$(campaign_name) — Inverse Scale (Flux Proxy): 1/L(z,t)",
        colorbar_title = " 1/L [m⁻¹]", opts...)

    # Panel 3: Dimensionless Stability Parameter ζ(z,t) = z/L
    p3 = heatmap(t_axis, z_axis, zeta_symlog;
        clims = (-1.5, 1.5), color = :pu_or,
        title = "$(campaign_name) — Stability Profile: sgn(ζ) log₁₀(1 + |ζ|)",
        colorbar_title = " sgn(ζ) log₁₀(1+|ζ|)", opts...)

    # Panel 4: Vertical Divergence ∂(1/L)/∂z
    d_invL_dz = zeros(Float64, nz, nt)
    for j in 1:nt
        d_invL_dz[:, j] = non_uniform_gradient_1d(data.inv_L[:, j], z_axis)
    end

    p4 = heatmap(t_axis, z_axis, clamp.(d_invL_dz, -0.01, 0.01);
        clims = (-0.005, 0.005), color = :vik,
        title = "$(campaign_name) — Vertical Scale Divergence: ∂(1/L)/∂z",
        colorbar_title = " ∂(1/L)/∂z [m⁻²]", opts...)

    fig = plot(p1, p2, p3, p4, layout = (2, 2), size = (1100, 750), dpi = 300)
    savefig(fig, out_path)
    return fig
end

# Helper non-uniform 1D derivative
function non_uniform_gradient_1d(f::Vector{Float64}, z::Vector{Float64})
    n = length(z)
    df = fill(0.0, n)
    length(z) < 3 && return df
    df[1] = (f[2] - f[1]) / (z[2] - z[1])
    df[end] = (f[end] - f[end-1]) / (z[end] - z[end-1])
    for i in 2:(n - 1)
        h1, h2 = z[i] - z[i-1], z[i+1] - z[i]
        df[i] = (f[i+1]*h1^2 - f[i-1]*h2^2 + f[i]*(h2^2 - h1^2)) / (h1 * h2 * (h1 + h2))
    end
    return df
end

# -------------------------------------------------------------------
# 3. Execution Pipeline Across Observational Datasets
# -------------------------------------------------------------------
function generate_all_sbltookit_heatmaps()
    for (name, path) in CAMPAIGN_SOURCES
        isfile(path) || continue
        df = CSV.read(path, DataFrame)

        # Ingest and pivot Obukhov fields
        obukhov_data = extract_obukhov_matrix(df)

        # Plot multi-panel observational comparison
        fig_out = joinpath(SBL_OUTPUT_DIR, "$(slugify(name))_obukhov_heatmaps.png")
        plot_sbltookit_obukhov_panel(obukhov_data, name, fig_out)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_all_sbltookit_heatmaps()
end