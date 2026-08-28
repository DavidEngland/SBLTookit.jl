#!/usr/bin/env julia
# src/visualize_ri_heatmaps.jl
# AUTHOR: David England (Refactored for Physical & GSPT Rigor)
# SBLToolkit Diagnostics & Visualization Suite
# ==========================================================================

using Pkg
Pkg.activate(dirname(@__DIR__))

using Dates
using Statistics

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

# --- 1. Smooth Mathematical Utilities ---

"""
    softplus_floor(x, x_min; k=10.0)

Smooth, infinitely differentiable approximation to max(x, x_min) using a log-sum-exp
formulation. Avoids C1 derivative discontinuities in finite-difference diagnostic chains.
"""
@inline function softplus_floor(x::T, x_min::T; k::T=T(10.0)) where {T<:AbstractFloat}
    return x_min + log1p(exp(k * (x - x_min))) / k
end


# --- 2. Synthetic SBL Diagnostic Generator ---

"""
    generate_synthetic_sbl_data()

Generates a synthetic 15-level idealized sparse SBL grid (z = 2 to 200m).
Uses differentiable profile functions to prevent artificial shear jump artifacts.
"""
function generate_synthetic_sbl_data()
    # 15-level idealized sparse SBL grid (Δz ≈ 14.1 m)
    z = collect(range(2.0, 200.0, length=15))
    nz = length(z)

    times = collect(range(0.0, 12.0, length=73)) # Hours since sunset
    nt = length(times)

    θ = zeros(nt, nz)
    u = zeros(nt, nz)
    v = zeros(nt, nz)

    θ0 = 265.0       # Reference surface temp (K)
    g = 9.81         # Acceleration due to gravity (m/s^2)
    h_θ = 45.0       # Inversion depth scale (m)
    w_jet = 35.0     # LLJ width scale (m)
    u_g = 4.5        # Geostrophic background wind (m/s)
    z0 = 0.1         # Roughness length (m)

    for i in 1:nt
        t = times[i]
        Δθ_surf = -12.0 * (1.0 - exp(-t / 3.5))

        for k in 1:nz
            # 1. Surface inversion profile
            θ[i, k] = θ0 + Δθ_surf * exp(-z[k] / h_θ) + 0.004 * z[k]

            # 2. Descending Low-Level Jet (LLJ)
            U_jet = 7.0 + 5.0 * sin(pi * t / 12.0)
            z_jet = 140.0 - 50.0 * (t / 12.0)

            u_jet_comp = U_jet * exp(-((z[k] - z_jet) / w_jet)^2)
            u_log_comp = u_g * (log(z[k] / z0) / log(200.0 / z0))

            u_raw = u_jet_comp + u_log_comp
            u[i, k] = softplus_floor(u_raw, 0.1; k=10.0)
            v[i, k] = 0.4 * u[i, k] * sin(pi * z[k] / 200.0)
        end
    end

    return times, z, θ, u, v, g
end


# --- 3. Robust Data Ingestion & Sentinel Normalization ---

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
    extract_netcdf_campaign(nc_file)

Flexible NetCDF reader with multi-campaign support (CASES-99, GABLS3, SCMs).
Dynamically infers spatial/temporal dimensions, scales units, and aligns matrices to (Time × Height).
"""
function extract_netcdf_campaign(nc_file::String)
    NCDataset(nc_file, "r") do ds
        # Helper for case-insensitive variable lookup
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
            error("Required profile variables (wind/temperature) not found in NetCDF schema.")
        end

        # Read native profile data arrays
        u_raw = Array(u_ds)
        θ_raw = Array(θ_ds)
        v_raw = v_ds !== nothing ? Array(v_ds) : zeros(Float64, size(u_raw))

        # Squeeze singleton dimensions if ndims > 2
        u_raw = dropdims(u_raw; dims=Tuple(i for i in 1:ndims(u_raw) if size(u_raw, i) == 1))
        θ_raw = dropdims(θ_raw; dims=Tuple(i for i in 1:ndims(θ_raw) if size(θ_raw, i) == 1))
        v_raw = dropdims(v_raw; dims=Tuple(i for i in 1:ndims(v_raw) if size(v_raw, i) == 1))

        # 1. Parse Time Axis
        if time_ds !== nothing
            t_data = Array(time_ds)
            if ndims(t_data) > 1
                t_data = t_data[:]
            end
            if eltype(t_data) <: Dates.AbstractTime
                times = (Dates.datetime2unix.(t_data) .- Dates.datetime2unix(t_data[1])) ./ 3600.0
            else
                times_raw = normalize_campaign_array(t_data)
                valid_t = filter(!isnan, times_raw)
                if !isempty(valid_t) && maximum(valid_t) > 1000.0
                    times_raw = (times_raw .- valid_t[1]) ./ 3600.0
                end
                times = times_raw
            end
        else
            times = Float64.(collect(1:size(u_raw, 1)))
        end
        nt = length(times)

        # 2. Infer nz & Align profile matrices to (nt × nz)
        s1, s2 = size(u_raw)
        if s1 == nt && s2 != nt
            nz = s2
            u_2d = normalize_campaign_array(u_raw)
            v_2d = normalize_campaign_array(v_raw)
            θ_2d = normalize_campaign_array(θ_raw)
        elseif s2 == nt && s1 != nt
            nz = s1
            u_2d = permutedims(normalize_campaign_array(u_raw), (2, 1))
            v_2d = permutedims(normalize_campaign_array(v_raw), (2, 1))
            θ_2d = permutedims(normalize_campaign_array(θ_raw), (2, 1))
        else
            # Fallback for square grids
            nz = (s1 == nt) ? s2 : s1
            u_2d = (s1 == nt) ? normalize_campaign_array(u_raw) : permutedims(normalize_campaign_array(u_raw), (2, 1))
            v_2d = (s1 == nt) ? normalize_campaign_array(v_raw) : permutedims(normalize_campaign_array(v_raw), (2, 1))
            θ_2d = (s1 == nt) ? normalize_campaign_array(θ_raw) : permutedims(normalize_campaign_array(θ_raw), (2, 1))
        end

        # 3. Resolve Vertical Height Coordinate (z) matching nz
        z_candidates = ["height", "Height", "level", "z", "altitude", "lev", "obs_height"]
        z = Float64[]

        for z_name in z_candidates
            var_cand = fetch_var([z_name])
            if var_cand !== nothing
                z_arr = Array(var_cand)
                # Handle 2D height arrays (nz × nt) or (nt × nz)
                if ndims(z_arr) == 2
                    if size(z_arr, 1) == nz
                        z_arr = vec(mean(z_arr, dims=2))
                    elseif size(z_arr, 2) == nz
                        z_arr = vec(mean(z_arr, dims=1))
                    end
                elseif ndims(z_arr) > 2
                    z_arr = vec(z_arr)
                end

                z_clean = filter(!isnan, normalize_campaign_array(z_arr))
                if length(z_clean) == nz
                    z = z_clean
                    break
                end
            end
        end

        # Fallback if no matching z variable of length nz was found
        if isempty(z)
            z = Float64.(collect(1:nz))
        end

        # 4. Temperature Unit Scaling (Celsius to Kelvin)
        valid_θ = filter(!isnan, θ_2d)
        if !isempty(valid_θ) && mean(valid_θ) < 100.0
            θ_2d = θ_2d .+ 273.15
        end

        return times, z, θ_2d, u_2d, v_2d, 9.81
    end
end


# --- 4. Diagnostic Calculation (Profiles & Derivatives) ---

"""
    process_campaign_profiles(times, z, θ, u, v, g; Ri_bound=2.0, S2_min=1e-8)

Computes buoyancy frequency N², shear S², raw Ri_g_raw, bounded Ri_g_reg,
regularization error E_Ri, vertical gradient ∂Ri/∂z, and curvature ∂²Ri/∂z².
"""
function process_campaign_profiles(times, z, θ, u, v, g; Ri_bound::Float64=2.0, S2_min::Float64=1e-8)
    nt = length(times)
    nz = length(z)

    z_mid = 0.5 .* (z[1:(end-1)] .+ z[2:end])
    nz_mid = length(z_mid)

    N2 = fill(NaN, nt, nz_mid)
    S2 = fill(NaN, nt, nz_mid)
    Ri_raw = fill(NaN, nt, nz_mid)
    Ri_reg = fill(NaN, nt, nz_mid)
    E_Ri = fill(NaN, nt, nz_mid)

    for i in 1:nt
        for k in 1:nz_mid
            Δz = z[k+1] - z[k]
            θ_mid = 0.5 * (θ[i, k] + θ[i, k+1])

            if !isnan(θ_mid) && θ_mid > 0 && !isnan(u[i, k]) && !isnan(u[i, k+1])
                θ_z = (θ[i, k+1] - θ[i, k]) / Δz
                u_z = (u[i, k+1] - u[i, k]) / Δz
                v_z = (v[i, k+1] - v[i, k]) / Δz

                N2[i, k] = (g / θ_mid) * θ_z
                S2[i, k] = u_z^2 + v_z^2

                if S2[i, k] > 0
                    Ri_raw[i, k] = N2[i, k] / S2[i, k]
                end

                S2_clamped = max(S2[i, k], S2_min)
                Ri_reg[i, k] = Ri_bound * tanh(N2[i, k] / (Ri_bound * S2_clamped))
                E_Ri[i, k] = Ri_reg[i, k] - Ri_raw[i, k]
            end
        end
    end

    # First derivative ∂Ri/∂z
    z_mid_mid = 0.5 .* (z_mid[1:(end-1)] .+ z_mid[2:end])
    Ri_z = fill(NaN, nt, nz_mid - 1)
    for i in 1:nt
        for k in 1:(nz_mid-1)
            Δz_mid = z_mid[k+1] - z_mid[k]
            if !isnan(Ri_reg[i, k+1]) && !isnan(Ri_reg[i, k])
                Ri_z[i, k] = (Ri_reg[i, k+1] - Ri_reg[i, k]) / Δz_mid
            end
        end
    end

    # Second derivative ∂²Ri/∂z² (Non-uniform central difference)
    z_curv = z_mid[2:(end-1)]
    Ri_zz = fill(NaN, nt, length(z_curv))
    for i in 1:nt
        for k in 2:(nz_mid-1)
            h1 = z_mid[k] - z_mid[k-1]
            h2 = z_mid[k+1] - z_mid[k]

            ri_prev = Ri_reg[i, k-1]
            ri_curr = Ri_reg[i, k]
            ri_next = Ri_reg[i, k+1]

            if !isnan(ri_prev) && !isnan(ri_curr) && !isnan(ri_next)
                dRi_dn = (ri_curr - ri_prev) / h1
                dRi_up = (ri_next - ri_curr) / h2
                Ri_zz[i, k-1] = 2.0 * (dRi_up - dRi_dn) / (h1 + h2)
            end
        end
    end

    return z_mid, z_mid_mid, z_curv, N2, S2, Ri_raw, Ri_reg, E_Ri, Ri_z, Ri_zz
end


# --- 5. Diagnostic Suite Visualization ---

"""
    create_gspt_diagnostic_figure(...)

Generates the 4-panel diagnostic figure.
"""
function create_gspt_diagnostic_figure(times, z_mid, z_mid_mid, z_curv, Ri_raw, Ri_reg, E_Ri, Ri_z, Ri_zz;
    campaign_name="15-level idealized sparse SBL grid",
    output_file="ri_gspt_diagnostic.png")

    fig = Figure(size=(1250, 1100), font="DejaVu Sans")
    Label(fig[0, 1:2], "Atmospheric Boundary Layer: Gradient Richardson Number Topology & Fold Dynamics",
        fontsize=16, font=:bold, halign=:left)

    # Panel A: Ri_reg
    ax1 = Axis(fig[1, 1], title="A. Regularized Richardson Field (Ri_g^reg)", ylabel="Height z (m)")
    hm1 = heatmap!(ax1, times, z_mid, Ri_reg, colormap=Reverse(:RdBu), colorrange=(-1.5, 1.5), nan_color=:gray90)
    try
        contour!(ax1, times, z_mid, Ri_reg, levels=[0.25], color=:black, linewidth=1.5, linestyle=:dash)
    catch
        ;
    end
    Colorbar(fig[1, 2], hm1, label="Ri_g^reg", ticks=-1.5:0.5:1.5)

    # Panel B: E_Ri Distortion
    ax2 = Axis(fig[2, 1], title="B. Visualization Regularization Bias (E_Ri = Ri_g^reg - Ri_g^raw)", ylabel="Height z (m)")
    hm2 = heatmap!(ax2, times, z_mid, E_Ri, colormap=:PuOr, colorrange=(-0.5, 0.5), nan_color=:gray90)
    Colorbar(fig[2, 2], hm2, label="ΔRi", ticks=-0.5:0.25:0.5)

    # Panel C: First Derivative ∂Ri/∂z
    ax3 = Axis(fig[3, 1], title="C. First Derivative (∂Ri_g^reg / ∂z) — Fold Line Tracking", ylabel="Height z (m)")
    hm3 = heatmap!(ax3, times, z_mid_mid, Ri_z, colormap=:curl, colorrange=(-0.05, 0.05), nan_color=:gray90)
    try
        contour!(ax3, times, z_mid_mid, Ri_z, levels=[0.0], color=:magenta, linewidth=2.0)
    catch
        ;
    end
    Colorbar(fig[3, 2], hm3, label="∂Ri / ∂z (m⁻¹)", ticks=-0.05:0.025:0.05)

    # Panel D: Vertical Curvature ∂²Ri/∂z²
    ax4 = Axis(fig[4, 1], title="D. Vertical Curvature (∂²Ri_g^reg / ∂z²) — Cusp Catastrophe Loci",
        xlabel="Time Axis", ylabel="Height z (m)")
    hm4 = heatmap!(ax4, times, z_curv, Ri_zz, colormap=:PRGn, colorrange=(-0.005, 0.005), nan_color=:gray90)
    try
        contour!(ax4, times, z_mid_mid, Ri_z, levels=[0.0], color=:magenta, linewidth=1.5)
        contour!(ax4, times, z_curv, Ri_zz, levels=[0.0], color=:orange, linewidth=1.5, linestyle=:dash)
    catch
        ;
    end
    Colorbar(fig[4, 2], hm4, label="∂²Ri / ∂z² (m⁻²)", ticks=-0.005:0.0025:0.005)

    Label(fig[5, 1:2], "Dataset: $campaign_name | Magenta: ∂Ri/∂z=0 (Folds), Orange: ∂²Ri/∂z²=0 (Inflections)",
        fontsize=10, color=:gray40, halign=:left, padding=(0, 0, 10, 0))

    save(output_file, fig, px_per_unit=2.0)
    println("SUCCESS: Graphic rendered to '$(output_file)'")
    return fig
end


# --- 6. Pipeline Execution ---

function main(args::Vector{String}=ARGS)
    println("==========================================================================")
    println("  GSPT SBL Diagnostic Suite: Richardson Number Topology & Fold Tracking  ")
    println("==========================================================================")

    if length(args) >= 1 && isfile(args[1])
        nc_file = args[1]
        println("Ingesting Campaign Dataset: '$nc_file'")
        times, z, θ, u, v, g = extract_netcdf_campaign(nc_file)
        camp_name = basename(nc_file)
        out_name = replace(camp_name, ".nc" => "_gspt_diagnostic.png")
    else
        println("Executing on 15-level idealized sparse SBL grid...")
        times, z, θ, u, v, g = generate_synthetic_sbl_data()
        camp_name = "15-level idealized sparse SBL grid"
        out_name = "ri_gspt_diagnostic.png"
    end

    z_mid, z_mid_mid, z_curv, N2, S2, Ri_raw, Ri_reg, E_Ri, Ri_z, Ri_zz = process_campaign_profiles(times, z, θ, u, v, g)

    create_gspt_diagnostic_figure(times, z_mid, z_mid_mid, z_curv, Ri_raw, Ri_reg, E_Ri, Ri_z, Ri_zz;
        campaign_name=camp_name, output_file=out_name)

    println("==========================================================================")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
else
    main(String[])
end