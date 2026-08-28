#!/usr/bin/env julia
# scripts/track_a_regularization.jl
# AUTHOR: David England (Refactored for Physical & GSPT Rigor)
# SBLToolkit Diagnostics: Track A (Primitive) vs Track B (Diagnostic) Regularization
# ==========================================================================

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dates
using Statistics
using LinearAlgebra
using Random
using Printf

try
    using NCDatasets
    using CairoMakie
catch e
    @error "Ensure dependencies are installed in project environment:"
    @info "julia --project=. -e 'using Pkg; Pkg.add([\"NCDatasets\", \"CairoMakie\"])'"
    rethrow(e)
end

using NCDatasets
using CairoMakie

# --- 1. Regularization Operators & Mathematical Utilities ---

"""
    build_discrete_roughness_matrix(z)

Constructs the discrete second-derivative roughness matrix (D2) for non-uniform vertical grids.
The resulting matrix R = D2' * D2 is a symmetric, positive semi-definite operator.
"""
function build_discrete_roughness_matrix(z::Vector{T}) where {T<:AbstractFloat}
    N = length(z)
    @assert N >= 3 "Grid must contain at least 3 levels to compute second derivatives."

    D2 = zeros(T, N - 2, N)
    for i in 1:(N-2)
        k = i + 1
        h1 = z[k] - z[k-1]
        h2 = z[k+1] - z[k]

        c_prev = T(2.0) / (h1 * (h1 + h2))
        c_curr = -T(2.0) / (h1 * h2)
        c_next = T(2.0) / (h2 * (h1 + h2))

        D2[i, k-1] = c_prev
        D2[i, k] = c_curr
        D2[i, k+1] = c_next
    end

    return D2' * D2
end

"""
    tikhonov_regularize_2d(M_raw, R, λ)

Applies 1D vertical Tikhonov regularization across every time step in a 2D matrix (nt × nz):
    (I + λ * R) * y_smooth = y_raw
"""
function tikhonov_regularize_2d(M_raw::Matrix{T}, R::Matrix{T}, λ::T) where {T<:AbstractFloat}
    nt, nz = size(M_raw)
    M_smooth = copy(M_raw)
    A = I + λ .* R

    for i in 1:nt
        row = M_raw[i, :]
        # Only solve if vector contains valid numeric data
        if !any(isnan, row)
            M_smooth[i, :] = A \ row
        else
            # Interpolate or fallback for missing data rows
            valid_idx = findall(!isnan, row)
            if length(valid_idx) >= 3
                row_filled = copy(row)
                # Quick nearest fill for missing boundaries before inversion
                for k in 1:nz
                    if isnan(row_filled[k])
                        closest = valid_idx[argmin(abs.(valid_idx .- k))]
                        row_filled[k] = row[closest]
                    end
                end
                M_smooth[i, :] = A \ row_filled
            end
        end
    end
    return M_smooth
end


# --- 2. Data Ingestion & Synthetic Campaign Generator ---

"""
    normalize_campaign_array(A; sentinels=[-9999.0, -999.0, 1e30])

Safely converts Julia `missing` values, NaNs, and sentinel flags into Float64 NaNs.
"""
function normalize_campaign_array(A::AbstractArray; sentinels=[-9999.0, -999.0, 1e30])
    A_clean = Array{Float64}(undef, size(A))
    for i in eachindex(A)
        val = A[i]
        if ismissing(val) || isnan(val) || any(s -> isapprox(Float64(val), s; rtol=1e-3), sentinels)
            A_clean[i] = NaN
        else
            A_clean[i] = Float64(val)
        end
    end
    return A_clean
end

"""
    generate_synthetic_llj_campaign(; σ_u=0.15, σ_θ=0.08, seed=42)

Generates a time-evolving 2D synthetic Low-Level Jet (LLJ) dataset with non-uniform tower spacing.
"""
function generate_synthetic_llj_campaign(; σ_u::Float64=0.15, σ_θ::Float64=0.08, seed::Int=42)
    Random.seed!(seed)
    z = Float64[2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 65.0, 80.0, 100.0, 120.0]
    nz = length(z)

    times = collect(range(0.0, 12.0, length=73)) # 12 Hours
    nt = length(times)

    u_obs = zeros(nt, nz)
    v_obs = zeros(nt, nz)
    θ_obs = zeros(nt, nz)

    θ0 = 265.0
    g = 9.81

    for i in 1:nt
        t = times[i]
        z_jet = 40.0 - 10.0 * sin(pi * t / 12.0) # Descending/oscillating jet nose
        u_max = 10.0 + 4.0 * sin(pi * t / 12.0)

        for k in 1:nz
            u_true = u_max * (z[k] / z_jet) * exp(1.0 - z[k] / z_jet)
            v_true = 0.2 * u_true
            θ_true = θ0 + 8.0 * tanh(z[k] / 50.0)

            u_obs[i, k] = u_true + σ_u * randn()
            v_obs[i, k] = v_true + σ_u * randn()
            θ_obs[i, k] = θ_true + σ_θ * randn()
        end
    end

    return times, z, θ_obs, u_obs, v_obs, g
end

"""
    extract_netcdf_campaign(nc_file)

Flexible NetCDF reader with campaign support (CASES-99, GABLS3, SCMs).
"""
function extract_netcdf_campaign(nc_file::String)
    NCDataset(nc_file, "r") do ds
        function fetch_var(names)
            for n in names
                if haskey(ds, n)
                    return ds[n]
                end
                for k in keys(ds)
                    if lowercase(k) == lowercase(n)
                        return ds[k]
                    end
                end
            end
            return nothing
        end

        time_ds = fetch_var(["time", "Time", "datetime", "base_time"])
        u_ds = fetch_var(["u", "U", "u_wind", "eastward_wind", "spd", "Spd"])
        v_ds = fetch_var(["v", "V", "v_wind", "northward_wind", "dir", "Dir"])
        θ_ds = fetch_var(["theta", "potential_temp", "pot_temp", "THETA", "temp", "T", "tdry", "tc"])

        if u_ds === nothing || θ_ds === nothing
            error("Required variables (wind/temperature) not found in NetCDF schema.")
        end

        u_raw = dropdims(Array(u_ds); dims=Tuple(i for i in 1:ndims(u_ds) if size(u_ds, i) == 1))
        θ_raw = dropdims(Array(θ_ds); dims=Tuple(i for i in 1:ndims(θ_ds) if size(θ_ds, i) == 1))
        v_raw = v_ds !== nothing ? dropdims(Array(v_ds); dims=Tuple(i for i in 1:ndims(v_ds) if size(v_ds, i) == 1)) : zeros(Float64, size(u_raw))

        if time_ds !== nothing
            t_data = vec(Array(time_ds))
            if eltype(t_data) <: Dates.AbstractTime
                times = (Dates.datetime2unix.(t_data) .- Dates.datetime2unix(t_data[1])) ./ 3600.0
            else
                times_raw = normalize_campaign_array(t_data)
                valid_t = filter(!isnan, times_raw)
                times = (!isempty(valid_t) && maximum(valid_t) > 1000.0) ? (times_raw .- valid_t[1]) ./ 3600.0 : times_raw
            end
        else
            times = Float64.(collect(1:size(u_raw, 1)))
        end
        nt = length(times)

        s1, s2 = size(u_raw)
        if s1 == nt && s2 != nt
            nz = s2
            u_2d, v_2d, θ_2d = normalize_campaign_array(u_raw), normalize_campaign_array(v_raw), normalize_campaign_array(θ_raw)
        else
            nz = s1
            u_2d = permutedims(normalize_campaign_array(u_raw), (2, 1))
            v_2d = permutedims(normalize_campaign_array(v_raw), (2, 1))
            θ_2d = permutedims(normalize_campaign_array(θ_raw), (2, 1))
        end

        z_candidates = ["height", "Height", "level", "z", "altitude", "lev", "obs_height"]
        z = Float64[]
        for z_name in z_candidates
            var_cand = fetch_var([z_name])
            if var_cand !== nothing
                z_arr = Array(var_cand)
                if ndims(z_arr) == 2
                    z_arr = size(z_arr, 1) == nz ? vec(mean(z_arr, dims=2)) : vec(mean(z_arr, dims=1))
                end
                z_clean = filter(!isnan, normalize_campaign_array(vec(z_arr)))
                if length(z_clean) == nz
                    z = z_clean
                    break
                end
            end
        end
        if isempty(z)
            z = Float64.(collect(1:nz))
        end

        valid_θ = filter(!isnan, θ_2d)
        if !isempty(valid_θ) && mean(valid_θ) < 100.0
            θ_2d .+= 273.15
        end

        return times, z, θ_2d, u_2d, v_2d, 9.81
    end
end


# --- 3. Track A vs Track B Pipeline Processing ---

"""
    compute_staggered_rig(u, v, θ, z, g; S2_min=2e-4)

Computes gradient Richardson profiles on staggered vertical mid-points.
"""
function compute_staggered_rig(u::Matrix{Float64}, v::Matrix{Float64}, θ::Matrix{Float64}, z::Vector{Float64}, g::Float64; S2_min::Float64=2e-4)
    nt, nz = size(u)
    z_mid = 0.5 .* (z[1:(end-1)] .+ z[2:end])
    nz_mid = length(z_mid)

    Ri = fill(NaN, nt, nz_mid)
    S2 = fill(NaN, nt, nz_mid)

    for i in 1:nt
        for k in 1:nz_mid
            Δz = z[k+1] - z[k]
            θ_mid = 0.5 * (θ[i, k] + θ[i, k+1])

            if !isnan(θ_mid) && θ_mid > 0 && !isnan(u[i, k]) && !isnan(u[i, k+1])
                u_z = (u[i, k+1] - u[i, k]) / Δz
                v_z = (v[i, k+1] - v[i, k]) / Δz
                θ_z = (θ[i, k+1] - θ[i, k]) / Δz

                N2 = (g / θ_mid) * θ_z
                shear2 = u_z^2 + v_z^2
                S2[i, k] = shear2

                if shear2 >= S2_min
                    Ri[i, k] = N2 / shear2
                end
                # S2 < S2_min remains NaN to isolate singularity
            end
        end
    end
    return z_mid, Ri, S2
end

"""
    process_regularization_tracks(times, z, θ_obs, u_obs, v_obs, g; λ_u=5.0, λ_θ=10.0, λ_Ri=5.0)

Executes primitive Track A vs diagnostic Track B comparison pipeline.
"""
function process_regularization_tracks(times, z, θ_obs, u_obs, v_obs, g; λ_u::Float64=5.0, λ_θ::Float64=10.0, λ_Ri::Float64=5.0)
    R_grid = build_discrete_roughness_matrix(z)

    # --- TRACK A: Regularize Primitive Fields First ---
    u_smooth_A = tikhonov_regularize_2d(u_obs, R_grid, λ_u)
    v_smooth_A = tikhonov_regularize_2d(v_obs, R_grid, λ_u)
    θ_smooth_A = tikhonov_regularize_2d(θ_obs, R_grid, λ_θ)

    z_mid, Ri_Track_A, S2_Track_A = compute_staggered_rig(u_smooth_A, v_smooth_A, θ_smooth_A, z, g; S2_min=2e-4)

    # --- TRACK B: Compute Raw Ri, then Regularize Diagnostic Field ---
    _, Ri_raw_B, _ = compute_staggered_rig(u_obs, v_obs, θ_obs, z, g; S2_min=2e-4)

    # Fill extreme NaNs/singularity spikes in raw diagnostic field prior to inversion
    Ri_raw_filled = copy(Ri_raw_B)
    for i in eachindex(Ri_raw_filled)
        if isnan(Ri_raw_filled[i]) || isinf(Ri_raw_filled[i])
            Ri_raw_filled[i] = 5.0 # Boundary fallback value
        end
    end

    R_mid = build_discrete_roughness_matrix(z_mid)
    Ri_Track_B = tikhonov_regularize_2d(Ri_raw_filled, R_mid, λ_Ri)

    # Quantify Singularity Leakage Bias: ΔRi = Track B - Track A
    Leakage_Bias = Ri_Track_B .- replace(Ri_Track_A, NaN => 0.0)

    return z_mid, Ri_Track_A, Ri_Track_B, Leakage_Bias, S2_Track_A
end


# --- 4. Diagnostic Visualization Suite ---

"""
    create_track_comparison_figure(...)

Generates 4-panel CairoMakie diagnostic plot comparing Track A and Track B.
"""
function create_track_comparison_figure(times, z_mid, Ri_Track_A, Ri_Track_B, Leakage_Bias, S2_Track_A;
    campaign_name="Idealized LLJ Campaign",
    output_file="track_a_regularization.png")

    fig = Figure(size=(1250, 1100), font="DejaVu Sans")
    Label(fig[0, 1:2], "Atmospheric Boundary Layer: Track A (Primitive) vs Track B (Diagnostic) Regularization",
        fontsize=16, font=:bold, halign=:left)

    # Panel A: Track A (Primitive Filtering)
    ax1 = Axis(fig[1, 1], title="A. Track A: Primitive Field Regularization (Ri_g from Smoothed u, v, θ)", ylabel="Height z (m)")
    hm1 = heatmap!(ax1, times, z_mid, Ri_Track_A, colormap=Reverse(:RdBu), colorrange=(-1.0, 2.0), nan_color=:gray90)
    try
        contour!(ax1, times, z_mid, Ri_Track_A, levels=[0.25], color=:black, linewidth=1.5, linestyle=:dash)
    catch
        ;
    end
    Colorbar(fig[1, 2], hm1, label="Ri_g (Track A)", ticks=-1.0:0.5:2.0)

    # Panel B: Track B (Direct Diagnostic Filtering)
    ax2 = Axis(fig[2, 1], title="B. Track B: Diagnostic Field Regularization (Direct Smoothing on Raw Ri_g)", ylabel="Height z (m)")
    hm2 = heatmap!(ax2, times, z_mid, Ri_Track_B, colormap=Reverse(:RdBu), colorrange=(-1.0, 2.0), nan_color=:gray90)
    Colorbar(fig[2, 2], hm2, label="Ri_g (Track B)", ticks=-1.0:0.5:2.0)

    # Panel C: Singularity Leakage Bias
    ax3 = Axis(fig[3, 1], title="C. Singularity Leakage Distortion Bias (ΔRi = Track B - Track A)", ylabel="Height z (m)")
    hm3 = heatmap!(ax3, times, z_mid, Leakage_Bias, colormap=:PuOr, colorrange=(-1.5, 1.5), nan_color=:gray90)
    Colorbar(fig[3, 2], hm3, label="ΔRi", ticks=-1.5:0.75:1.5)

    # Panel D: Shear Field & Jet Core Loci
    ax4 = Axis(fig[4, 1], title="D. Primitive Shear Field Squared (S²) — Low-Level Jet Core Singularity Region",
        xlabel="Time Axis", ylabel="Height z (m)")
    hm4 = heatmap!(ax4, times, z_mid, S2_Track_A, colormap=:magma, colorrange=(0.0, 0.01), nan_color=:gray90)
    try
        contour!(ax4, times, z_mid, S2_Track_A, levels=[2e-4], color=:cyan, linewidth=1.5, linestyle=:dash)
    catch
        ;
    end
    Colorbar(fig[4, 2], hm4, label="S² (s⁻²)", ticks=0.0:0.0025:0.01)

    Label(fig[5, 1:2], "Dataset: $campaign_name | Cyan Dash: S² Minimum Floor (Singularity Boundary)",
        fontsize=10, color=:gray40, halign=:left, padding=(0, 0, 10, 0))

    save(output_file, fig, px_per_unit=2.0)
    println("SUCCESS: Graphic rendered to '$(output_file)'")
    return fig
end


# --- 5. Pipeline Execution ---

function main(args::Vector{String}=ARGS)
    println("==========================================================================")
    println("  GSPT SBL Diagnostics: Track A (Primitive) Regularization Pipeline      ")
    println("==========================================================================")

    if length(args) >= 1 && isfile(args[1])
        nc_file = args[1]
        println("Ingesting Campaign Dataset: '$nc_file'")
        times, z, θ_obs, u_obs, v_obs, g = extract_netcdf_campaign(nc_file)
        camp_name = basename(nc_file)
        out_name = replace(camp_name, ".nc" => "_track_a_comparison.png")
    else
        println("Executing on Synthetic Low-Level Jet (LLJ) Campaign...")
        times, z, θ_obs, u_obs, v_obs, g = generate_synthetic_llj_campaign()
        camp_name = "Synthetic Low-Level Jet Campaign"
        out_name = "track_a_regularization.png"
    end

    z_mid, Ri_Track_A, Ri_Track_B, Leakage_Bias, S2_Track_A = process_regularization_tracks(times, z, θ_obs, u_obs, v_obs, g)

    create_track_comparison_figure(times, z_mid, Ri_Track_A, Ri_Track_B, Leakage_Bias, S2_Track_A;
        campaign_name=camp_name, output_file=out_name)

    println("==========================================================================")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
else
    main(String[])
end