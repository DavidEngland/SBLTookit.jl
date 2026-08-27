#!/usr/bin/env julia
# src/visualize_ri_heatmaps.jl
# AUTHOR: David England
# SBLToolkit Visualization Script
# ==========================================================================
# FILE: visualize_ri_heatmaps.jl
# DESCRIPTION: High-fidelity visualization script using Julia (CairoMakie + NCDatasets)
#              to compute and map regularized Richardson number (Ri_g) profiles
#              as spatial-temporal heatmaps across atmospheric campaign datasets.
#              Supports both standard NetCDF campaign ingestion and a robust,
#              physically-derived Stable Boundary Layer (SBL) synthetic demo.
# ==========================================================================

using Dates
using Statistics

# --- 1. Environment & Package Check ---
# We use standard Julia package loading blocks with clear warnings.
# Since Julia packages might need to be pre-installed by the user, we handle
# imports defensively and provide instructions.
try
    using NCDatasets
    using CairoMakie
catch e
    @error "Required visualization packages (NCDatasets.jl, CairoMakie.jl) not fully loaded."
    @info "To run this script, please ensure they are installed in your Julia environment:"
    @info "julia> using Pkg; Pkg.add([\"NCDatasets\", \"CairoMakie\"])"
    rethrow(e)
end

"""
    compute_safe_rig(N2, S2; Ri_bound=2.0, ϵ_s=1e-12)

Computes the regularized, safe, and smoothly bounded gradient Richardson number (Ri_g).
Utilizes:
1. Denominator regularization floor (ϵ_s) to prevent division-by-zero under vanishing shear.
2. Asymptotic hyperbolic tangent clamp (Ri_bound) to elegantly bound extremely stable or
   unstable convective regimes, preventing solver stiffness and providing a visually clean range.
"""
@inline function compute_safe_rig(N2::T, S2::T; Ri_bound::T=T(2.0), ϵ_s::T=T(1e-12)) where {T <: AbstractFloat}
    # Regularize vertical shear to prevent singularity
    S2_safe = S2 + ϵ_s
    # Apply analytic hyperbolic tanh clamp to map raw range (-inf, +inf) smoothly to [-Ri_bound, +Ri_bound]
    return Ri_bound * tanh(N2 / (Ri_bound * S2_safe))
end

# Handle generic numbers by promotion to Float64
compute_safe_rig(N2::Real, S2::Real; Ri_bound::Real=2.0, ϵ_s::Real=1e-12) = 
    compute_safe_rig(Float64(N2), Float64(S2); Ri_bound=Float64(Ri_bound), ϵ_s=Float64(ϵ_s))


"""
    generate_synthetic_sbl_data()

Generates a physically-grounded nocturnal Stable Boundary Layer (SBL) scenario to serve as
a high-fidelity diagnostic testbed. Simulates:
1. Progressive surface radiative cooling after sunset, producing a steep potential temperature
   inversion near the ground (high θ_z) that weakens with height.
2. An evolving Low-Level Jet (LLJ) in the u-wind component that intensifies and slowly descends
   over the 12-hour nocturnal period, generating strong shear layers below the jet core.
3. Ekman-friction-induced wind turning in the v-wind component.
"""
function generate_synthetic_sbl_data()
    # 15 vertical levels (standardized campaign grid, e.g., SHEBA/FLOSS-like tower levels)
    z = collect(range(2.0, 200.0, length=15))
    nz = length(z)
    
    # 12-hour night with 10-minute sampling resolution (73 time points)
    times = collect(range(0.0, 12.0, length=73)) # Hours since sunset
    nt = length(times)
    
    # Pre-allocate 2D profiles (Time x Height)
    θ = zeros(nt, nz)
    u = zeros(nt, nz)
    v = zeros(nt, nz)
    
    # Physical Constants
    θ0 = 265.0       # Cold winter reference temperature (K)
    g = 9.81         # Gravitational acceleration (m/s^2)
    h_θ = 45.0       # Inversion depth scale (m)
    w_jet = 35.0     # Low-level jet width scale (m)
    u_g = 4.5        # Background geostrophic wind (m/s)
    
    for i in 1:nt
        t = times[i]
        
        # A. Radiative Cooling: Surface potential temperature cools asymptotically with time
        Δθ_surf = -12.0 * (1.0 - exp(-t / 3.5)) # Up to 12 K cooling at surface
        for k in 1:nz
            # Exponential inversion profile + background weak tropospheric stability (0.004 K/m)
            θ[i, k] = θ0 + Δθ_surf * exp(-z[k] / h_θ) + 0.004 * z[k]
        end
        
        # B. Descending Low-Level Jet (LLJ): Peak wind speed intensifies at midnight, core descends
        U_jet = 7.0 + 5.0 * sin(pi * t / 12.0)        # Jet speed ranges from 7m/s to 12m/s
        z_jet = 140.0 - 50.0 * (t / 12.0)             # Jet core descends from 140m to 90m over night
        
        for k in 1:nz
            # Gaussian jet structure + logarithmic boundary layer recovery
            u_jet_comp = U_jet * exp(-((z[k] - z_jet) / w_jet)^2)
            u_log_comp = u_g * (log(z[k] / 0.1) / log(200.0 / 0.1))
            
            u[i, k] = max(0.1, u_jet_comp + u_log_comp)
            # Coriolis-induced directional wind turning in the SBL
            v[i, k] = 0.4 * u[i, k] * sin(pi * z[k] / 200.0)
        end
    end
    
    return times, z, θ, u, v, θ0, g
end


"""
    process_campaign_profiles(times, z, θ, u, v, θ0, g; Ri_bound=2.0, ϵ_s=1e-12)

Computes buoyancy gradients, horizontal wind shears, and safe/smooth gradient Richardson profiles
on the vertical mid-levels between tower height coordinates.
"""
function process_campaign_profiles(times, z, θ, u, v, θ0, g; Ri_bound=2.0, ϵ_s=1e-12)
    nt = length(times)
    nz = length(z)
    
    # Calculate mid-level heights for staggered vertical gradients
    z_mid = 0.5 .* (z[1:end-1] .+ z[2:end])
    nz_mid = length(z_mid)
    
    # Pre-allocate gradient arrays
    N2 = zeros(nt, nz_mid)
    S2 = zeros(nt, nz_mid)
    Ri_g = zeros(nt, nz_mid)
    
    g_over_θ0 = g / θ0
    
    for i in 1:nt
        for k in 1:nz_mid
            Δz = z[k+1] - z[k]
            
            # Finite-difference gradients on staggered grid levels
            θ_z = (θ[i, k+1] - θ[i, k]) / Δz
            u_z = (u[i, k+1] - u[i, k]) / Δz
            v_z = (v[i, k+1] - v[i, k]) / Δz
            
            # Buoyancy frequency squared (buoyancy gradient)
            N2[i, k] = g_over_θ0 * θ_z
            
            # Vertical wind shear squared
            S2[i, k] = u_z^2 + v_z^2
            
            # Compute safe & smoothly clamped Richardson number
            # If values are missing/NaN, populate with NaN
            if isnan(N2[i, k]) || isnan(S2[i, k])
                Ri_g[i, k] = NaN
            else
                Ri_g[i, k] = compute_safe_rig(N2[i, k], S2[i, k]; Ri_bound=Ri_bound, ϵ_s=ϵ_s)
            end
        end
    end
    
    return z_mid, N2, S2, Ri_g
end


"""
    create_ri_heatmap(times, z_mid, Ri_g; campaign_name="Synthetic SBL", output_file="ri_heatmap.png")

Renders a publication-quality spatial-temporal heatmap of the regularized Richardson number.
Employs an elegant diverging color scale, highlights the critical threshold (Ri_c = 0.25) with a 
dashed contour line, and structures labels to convey immediate physical insights.
"""
function create_ri_heatmap(times, z_mid, Ri_g; campaign_name="Synthetic SBL", output_file="ri_heatmap.png")
    # Setup Figure and custom thematic elements
    # CairoMakie guarantees crisp vector and high-DPI rasterization
    fig = Figure(size = (1000, 650), font = "DejaVu Sans")
    
    # Title-as-takeaway: Tells the scientific story
    title_str = "Nocturnal Stable Boundary Layer: Low-Level Jet Shear Suppresses Richardson Number below Jet Core"
    if campaign_name != "Synthetic SBL"
        title_str = "Campaign $(campaign_name): Spatial-Temporal Gradient Richardson Number Profile"
    end
    
    ax = Axis(fig[1, 1],
        title = title_str,
        xlabel = campaign_name == "Synthetic SBL" ? "Hours Since Sunset (t)" : "Time Index / Date",
        ylabel = "Height Above Ground z (m)",
        titlesize = 14,
        titlealign = :left,
        titlefont = :bold,
        xgridstyle = :dash, ygridstyle = :dash,
        xgridcolor = :gray90, ygridcolor = :gray90,
        xminorticks = IntervalsBetween(5),
        yminorticks = IntervalsBetween(5),
        xminorticksvisible = true, yminorticksvisible = true
    )
    
    # 1. Plot Heatmap
    # We map the bounded Richardson number on a diverging :RdBu colormap (reversed to make 
    # blue represent stable stratification Ri_g > 0 and red represent unstable mixing Ri_g < 0).
    # Since stable stratification (blue) and unstable shear (red) are physical opposites, 
    # the neutral threshold (0.0) aligns perfectly with the white/subdued midpoint.
    hm = heatmap!(ax, times, z_mid, Ri_g;
        colormap = Reverse(:RdBu),
        colorrange = (-1.5, 1.5), # Highlight detailed structure around neutral-to-moderately-stable regimes
        nan_color = :gray93
    )
    
    # 2. Add Critical Boundary Contour (Ri_c = 0.25)
    # The critical Richardson number represents the physical transition where the turbulent engine dies.
    # We superimpose a dashed black contour line to mark this exact boundary.
    try
        contour!(ax, times, z_mid, Ri_g;
            levels = [0.25],
            color = :black,
            linewidth = 2.0,
            linestyle = :dash
        )
        # Add visual legend helper in a text box
        text!(ax, minimum(times) + 0.5, maximum(z_mid) - 15.0, 
            text = "--- Ri_c = 0.25 (Turbulent Transition)", 
            color = :black, 
            fontsize = 11,
            font = :bold
        )
    catch
        @warn "Contour plotting skipped due to data bounds or package limitation."
    end
    
    # 3. Add Annotation Callouts for Key SBL Dynamics
    if campaign_name == "Synthetic SBL"
        # Label the unstable ground layer under surface warming or high shear
        text!(ax, 2.0, 10.0, text = "Shear-Generated Turbulence (Ri < 0.25)", color = :darkred, fontsize = 10, font = :bold)
        # Label the supercritical laminar layer above the jet where mixing ceases
        text!(ax, 8.0, 160.0, text = "Supercritical Stratified Layer (Ri > Ri_c)", color = :darkblue, fontsize = 10, font = :bold)
        # Label descending jet core shear trace
        text!(ax, 5.0, 85.0, text = "Desending LLJ Shear Zone", color = :black, fontsize = 10, font = :italic)
    end
    
    # 4. Add Colorbar with Clean Interval Ticks
    cb = Colorbar(fig[1, 2], hm,
        label = "Regularized Gradient Richardson Number (Ri_g)",
        labelsize = 12,
        ticks = (-1.5:0.5:1.5, ["<-1.5 (Convective)", "-1.0", "-0.5", "0.0 (Neutral)", "0.5 (Stable)", "1.0", ">1.5 (Laminar)"]),
        width = 20,
        ticklabelsize = 10
    )
    
    # 5. Add Source Footnote
    Label(fig[2, 1], "Source: Campaign Dataset Archive, Z0HR Safe & Smooth Regularization (ϵ_s = 1e-12, Ri_bound = 2.0)",
        fontsize = 10, color = :gray50, halign = :left, padding = (0, 0, 10, 0))
    
    # Apply tight layout padding and save high-fidelity output
    rowsize!(fig.layout, 1, Relative(0.92))
    save(output_file, fig, px_per_unit = 2.0) # Double pixel resolution (300 DPI equivalent)
    println("SUCCESS: High-fidelity spatial-temporal heatmap saved to '$(output_file)'")
    
    return fig
end


"""
    main(args)

Main execution routine. If a NetCDF filepath is passed in `args[1]`, the script parses 
the observational variables. Otherwise, it triggers the synthetic SBL testbed.
"""
function main(args::Vector{String} = ARGS)
    println("==========================================================================")
    println("   Z0HR Safe & Smooth Richardson Heatmap Ingestion & Visualization Suite  ")
    println("==========================================================================")
    
    # Define hyperparameter defaults
    Ri_bound = 2.0
    ϵ_s = 1e-12
    
    if length(args) >= 1 && isfile(args[1])
        nc_file = args[1]
        println("Ingesting Campaign NetCDF Dataset: '$(nc_file)'...")
        
        try
            NCDataset(nc_file, "r") do ds
                println("Checking NetCDF schema attributes...")
                # Fetch dimensions and coordinates
                times = haskey(ds, "time") ? ds["time"][:] : collect(1:size(ds["u"], 1))
                z = haskey(ds, "height") ? ds["height"][:] : (haskey(ds, "level") ? ds["level"][:] : collect(1:size(ds["u"], 2)))

                # NetCDF time dims are often DateTime; convert to elapsed hours for plotting/arithmetic
                if eltype(times) <: Dates.AbstractTime
                    times = Dates.datetime2unix.(times) .- Dates.datetime2unix(times[1])
                    times = times ./ 3600.0
                end
                
                # Fetch raw physical variables
                u = ds["u"][:, :]
                v = ds["v"][:, :]
                θ = haskey(ds, "theta") ? ds["theta"][:, :] : (haskey(ds, "potential_temp") ? ds["potential_temp"][:, :] : ds["temp"][:, :])
                
                # Deduce metadata reference fields
                g = 9.81
                θ0 = mean(filter(!isnan, θ))
                
                println("Computing gradients across $(length(times)) timestamps and $(length(z)) vertical levels...")
                z_mid, N2, S2, Ri_g = process_campaign_profiles(times, z, θ, u, v, θ0, g; Ri_bound=Ri_bound, ϵ_s=ϵ_s)
                
                output_name = replace(basename(nc_file), ".nc" => "_ri_heatmap.png")
                create_ri_heatmap(times, z_mid, Ri_g; campaign_name=basename(nc_file), output_file=output_name)
            end
        catch e
            @error "Failed to parse NetCDF campaign file. Error: $e"
            println("Falling back to SBL Synthetic testbed execution...")
            run_synthetic_flow(Ri_bound, ϵ_s)
        end
    else
        println("No valid NetCDF campaign file provided as argument.")
        println("Executing high-fidelity physical SBL synthetic demonstration...")
        run_synthetic_flow(Ri_bound, ϵ_s)
    end
    
    println("==========================================================================")
    return 0
end

function run_synthetic_flow(Ri_bound, ϵ_s)
    times, z, θ, u, v, θ0, g = generate_synthetic_sbl_data()
    println("Processing SBL gradients (15 standardized levels, 73 contiguous intervals)...")
    z_mid, N2, S2, Ri_g = process_campaign_profiles(times, z, θ, u, v, θ0, g; Ri_bound=Ri_bound, ϵ_s=ϵ_s)
    output_path = joinpath(pwd(), "ri_temporal_spatial_heatmap.png")
    create_ri_heatmap(times, z_mid, Ri_g; campaign_name="Synthetic SBL", output_file=output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
else
    # Automatically execute synthetic demo when evaluated in notebook/REPL context
    run_synthetic_flow(2.0, 1e-12)
end
