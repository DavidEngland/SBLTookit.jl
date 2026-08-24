#!/usr/bin/env julia
# SBLToolkit.jl/scripts/plot_obukhov_heatmaps.jl

using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

using SBLToolkit
using SBLToolkit.DatasetRegistry
using CSV
using DataFrames
using Dates
using NCDatasets
using Plots
using Plots.PlotMeasures
using Statistics

const SBL_OUTPUT_DIR = joinpath(PROJECT_ROOT, "reports", "generated", "sbltoolkit_heatmaps")
mkpath(SBL_OUTPUT_DIR)

slugify(s::String) = lowercase(replace(s, r"[^A-Za-z0-9]+" => "_"))
norm_str(s::Union{String,Symbol}) = lowercase(replace(string(s), r"[^a-z0-9]" => ""))

function safe_float(val)
    val === nothing && return NaN
    val === missing && return NaN
    val isa Number && return Float64(val)
    if val isa AbstractString
        parsed = tryparse(Float64, strip(val))
        return parsed !== nothing ? parsed : NaN
    end
    return NaN
end

# Convert local Gradient Richardson number Ri_g to stability parameter zeta
function rig_to_zeta(rig::Float64)
    !isfinite(rig) && return NaN
    if rig >= 0.0
        rig_c = min(rig, 0.19) # Cap below critical Ri_g ~ 0.2 singularity
        return rig_c / (1.0 - 5.0 * rig_c)
    else
        return rig
    end
end

function parse_rig_height(col_name::Symbol)
    s = string(col_name)
    m = match(r"^ri_g_(\d+)_(\d+)$"i, s)
    if m !== nothing
        i_part = parse(Float64, m[1])
        d_part = parse(Float64, m[2])
        return i_part + d_part / (10^length(m[2]))
    end
    m2 = match(r"^ri_g_(\d+)$"i, s)
    if m2 !== nothing
        return parse(Float64, m2[1])
    end
    return nothing
end

# -------------------------------------------------------------------
# Two-Pass Target Matching
# -------------------------------------------------------------------
function match_fuzzy_col(cols::Vector{Symbol}, targets::Vector{String})
    norm_targets = norm_str.(targets)
    norm_cols = norm_str.(cols)

    for target in norm_targets
        for (i, nc) in enumerate(norm_cols)
            nc == target && return cols[i]
        end
    end

    for target in norm_targets
        length(target) < 3 && continue
        for (i, nc) in enumerate(norm_cols)
            occursin(target, nc) && return cols[i]
        end
    end

    return nothing
end

function match_nc_var(keys_list::Vector{String}, targets::Vector{String})
    norm_targets = norm_str.(targets)
    norm_keys = norm_str.(keys_list)

    for target in norm_targets
        for (i, nk) in enumerate(norm_keys)
            nk == target && return keys_list[i]
        end
    end

    for target in norm_targets
        length(target) < 3 && continue
        for (i, nk) in enumerate(norm_keys)
            occursin(target, nk) && return keys_list[i]
        end
    end

    return nothing
end

# -------------------------------------------------------------------
# Stencil Differential Operators
# -------------------------------------------------------------------
function non_uniform_gradient_1d(f::Vector{Float64}, z::Vector{Float64})
    n = length(z)
    df = fill(NaN, n)
    valid_idx = findall(isfinite, f)
    length(valid_idx) < 3 && return df

    zv, fv = z[valid_idx], f[valid_idx]
    nv = length(zv)
    dfv = zeros(Float64, nv)

    dfv[1] = (fv[2] - fv[1]) / (zv[2] - zv[1])
    dfv[end] = (fv[end] - fv[end-1]) / (zv[end] - zv[end-1])

    for i in 2:(nv-1)
        h1, h2 = zv[i] - zv[i-1], zv[i+1] - zv[i]
        dfv[i] = (fv[i+1]*h1^2 - fv[i-1]*h2^2 + fv[i]*(h2^2 - h1^2)) / (h1 * h2 * (h1 + h2))
    end

    df[valid_idx] .= dfv
    return df
end

function non_uniform_hessian_1d(f::Vector{Float64}, z::Vector{Float64})
    return non_uniform_gradient_1d(non_uniform_gradient_1d(f, z), z)
end

function assemble_derivatives_from_zeta(timestamps, z_levels, zeta_mat)
    nz, nt = length(z_levels), length(timestamps)

    inv_L_mat = fill(NaN, nz, nt)
    for j in 1:nt, i in 1:nz
        if isfinite(zeta_mat[i, j]) && z_levels[i] > 1e-3
            inv_L_mat[i, j] = zeta_mat[i, j] / z_levels[i]
        end
    end

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
# Trajectory CSV Parsing Pipeline
# -------------------------------------------------------------------
function parse_obukhov_from_df(df::DataFrame, campaign_name::String)
    cols = propertynames(df)

    time_col = match_fuzzy_col(cols, ["sample_index", "sampleindex", "time_value", "time", "datetime", "index"])
    time_col === nothing && error("No valid time column found in trajectory dataframe.")

    timestamps = sort(unique([safe_float(v) for v in df[!, time_col]]))
    timestamps = filter(isfinite, timestamps)
    nt = length(timestamps)
    time_map = Dict(t => j for (j, t) in enumerate(timestamps))

    # Mode 1: Wide Gradient Richardson profiles (ri_g_X_Y)
    rig_cols = Symbol[]
    rig_heights = Float64[]
    for c in cols
        h = parse_rig_height(c)
        if h !== nothing
            push!(rig_cols, c)
            push!(rig_heights, h)
        end
    end

    if length(rig_cols) >= 3
        p = sortperm(rig_heights)
        z_levels = rig_heights[p]
        sorted_rig_cols = rig_cols[p]
        nz = length(z_levels)
        zeta_mat = fill(NaN, nz, nt)

        for (i, col) in enumerate(sorted_rig_cols)
            for row in eachrow(df)
                t_val = safe_float(row[time_col])
                if haskey(time_map, t_val)
                    rig_val = safe_float(row[col])
                    zeta_mat[i, time_map[t_val]] = rig_to_zeta(rig_val)
                end
            end
        end

        @info "Parsed Multi-Level Ri_g profiles [$campaign_name]: $(count(isfinite, zeta_mat))/$(length(zeta_mat)) non-NaN entries across $(nz) heights."
        return assemble_derivatives_from_zeta(timestamps, z_levels, zeta_mat)
    end

    # Mode 2: Standard 1D L_obukhov column fallback
    l_col = match_fuzzy_col(cols, ["L_obukhov", "lobukhov", "obukhovlength", "l"])
    if l_col !== nothing
        z_levels = [1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 120.0, 180.0]
        nz = length(z_levels)
        zeta_mat = fill(NaN, nz, nt)

        for row in eachrow(df)
            t_val = safe_float(row[time_col])
            if haskey(time_map, t_val)
                L_val = safe_float(row[l_col])
                if isfinite(L_val) && abs(L_val) > 1e-4
                    for i in 1:nz
                        zeta_mat[i, time_map[t_val]] = z_levels[i] / L_val
                    end
                end
            end
        end

        @info "Parsed 1D Surface L_obukhov [$campaign_name]: $(count(isfinite, zeta_mat))/$(length(zeta_mat)) non-NaN entries."
        return assemble_derivatives_from_zeta(timestamps, z_levels, zeta_mat)
    end

    error("No usable multi-level Ri_g profiles or L_obukhov columns in trajectory file.")
end

# -------------------------------------------------------------------
# NetCDF Loader Fallback
# -------------------------------------------------------------------
function search_nc_file(campaign_name::String)
    clean = replace(lowercase(campaign_name), r"[^a-z0-9]" => "")
    candidate_dirs = [
        joinpath(PROJECT_ROOT, "data", clean),
        joinpath(PROJECT_ROOT, "data", "raw", clean),
        joinpath(PROJECT_ROOT, "data")
    ]

    for d in candidate_dirs
        !isdir(d) && continue
        for (root, _, files) in walkdir(d)
            for f in files
                if endswith(f, ".nc") && (occursin("prof", lowercase(f)) || occursin("main", lowercase(f)) || occursin(clean, lowercase(f)))
                    return joinpath(root, f)
                end
            end
        end
    end
    return nothing
end

function load_profile_from_nc(nc_path::String, campaign_name::String; g=9.81, theta_0=273.15, k_vk=0.40, rho_cp=1200.0)
    !isfile(nc_path) && return nothing

    return NCDataset(nc_path, "r") do ds
        keys_list = collect(keys(ds))

        t_key = match_nc_var(keys_list, ["time", "datetime", "t", "sampleindex"])
        t_key === nothing && error("Missing time dimension in $nc_path.")
        timestamps = Float64.(vec(ds[t_key][:]))

        z_levels = [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 80.0, 100.0, 180.0]
        u_star_key = match_nc_var(keys_list, ["ustar", "ustarm", "frictionvelocity"])
        h_key = match_nc_var(keys_list, ["hs", "heatflux", "shf", "wt", "kinematicheatflux", "h", "qsurface"])

        (u_star_key === nothing || h_key === nothing) && error("Missing flux keys in NetCDF.")

        u_s = vec(Float64.(ds[u_star_key][:]))
        raw_h = vec(Float64.(ds[h_key][:]))
        nt, nz = length(timestamps), length(z_levels)
        zeta_mat = fill(NaN, nz, nt)

        for j in 1:nt
            wt = abs(raw_h[j]) > 5.0 ? raw_h[j] / rho_cp : raw_h[j]
            if isfinite(u_s[j]) && isfinite(wt) && abs(wt) > 1e-5
                L_val = -(u_s[j]^3 * theta_0) / (k_vk * g * wt)
                for i in 1:nz
                    zeta_mat[i, j] = z_levels[i] / L_val
                end
            end
        end

        @info "NetCDF Extracted Flux Profiles [$campaign_name]: $(count(isfinite, zeta_mat))/$(length(zeta_mat)) non-NaN entries."
        return assemble_derivatives_from_zeta(timestamps, z_levels, zeta_mat)
    end
end

# -------------------------------------------------------------------
# Plotting Routine
# -------------------------------------------------------------------
function plot_sbltoolkit_obukhov_panel(data::NamedTuple, campaign_name::String, out_path::String)
    t_axis = data.timestamps
    z_axis = data.z_levels

    opts = (ylabel="Height z [m]", framestyle=:box, guidefontsize=9, tickfontsize=8, titlefontsize=10, margin=3mm)

    cg_bwr = cgrad([:blue, :white, :red])
    cg_puor = cgrad([:purple, :white, :orange])

    p1 = heatmap(t_axis, z_axis, clamp.(data.inv_L, -0.1, 0.1);
        clims=(-0.05, 0.05), color=cg_bwr, title="$(campaign_name) — 1/L(z,t)", colorbar_title=" 1/L [m⁻¹]", opts...)

    zeta_symlog = @. sign(data.zeta) * log10(1.0 + abs(data.zeta))
    p2 = heatmap(t_axis, z_axis, zeta_symlog;
        clims=(-1.5, 1.5), color=cg_puor, title="$(campaign_name) — ζ(z,t)", colorbar_title=" sgn(ζ) log₁₀(1+|ζ|)", opts...)

    p3 = heatmap(t_axis, z_axis, clamp.(data.zeta_z, -0.2, 0.2);
        clims=(-0.1, 0.1), color=cg_bwr, title="$(campaign_name) — Jacobian ζ_z", colorbar_title=" ζ_z [m⁻¹]", opts...)

    if count(isfinite, data.zeta_z) > 0
        contour!(p3, t_axis, z_axis, data.zeta_z, levels=[0.0], color=:white, lw=1.5, ls=:solid)
    end

    p4 = heatmap(t_axis, z_axis, clamp.(data.zeta_zz, -0.05, 0.05);
        clims=(-0.02, 0.02), color=cg_bwr, title="$(campaign_name) — Curvature ζ_zz", colorbar_title=" ζ_zz [m⁻²]", opts...)

    fig = plot(p1, p2, p3, p4, layout=(2, 2), size=(1100, 750), dpi=300)
    savefig(fig, out_path)
    return fig
end

# -------------------------------------------------------------------
# Driver
# -------------------------------------------------------------------
function main()
    campaign_sources = discover_campaign_trajectories()
    @info "Discovered $(length(campaign_sources)) campaign trajectories from manifest schema."

    success_count = 0

    for (name, path) in campaign_sources
        @info "Processing $name..."
        obukhov_data = nothing

        # 1. Parse CSV trajectory file
        if isfile(path)
            try
                df = CSV.read(path, DataFrame)
                obukhov_data = parse_obukhov_from_df(df, name)
            catch e
                @warn "Trajectory CSV parsing failed for $name: $e"
            end
        end

        # 2. NetCDF fallback
        if obukhov_data === nothing
            nc_file = search_nc_file(name)
            if nc_file !== nothing
                @info "Attempting NetCDF extraction for $name via $nc_file"
                try
                    obukhov_data = load_profile_from_nc(nc_file, name)
                catch e
                    @warn "NetCDF fallback failed for $name: $e"
                end
            end
        end

        if obukhov_data === nothing
            @warn "Skipping $name: No valid stability profile data found."
            continue
        end

        fig_out = joinpath(SBL_OUTPUT_DIR, "$(slugify(name))_obukhov_heatmaps.png")
        plot_sbltoolkit_obukhov_panel(obukhov_data, name, fig_out)
        @info "Successfully rendered heatmap for $name -> $fig_out"
        success_count += 1
    end

    @info "Completed heatmaps pipeline: $success_count/$(length(campaign_sources)) rendered."
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end