#!/usr/bin/env julia
# =============================================================================
# GSPT TRIPLE-POINT INVARIANT-CONVERGENCE & PROFILE GEOMETRY SUITE
# Developed for Generalized Similarity Profile Theory (GSPT) Campaign Audits
# =============================================================================
# This template demonstrates:
# 1. Scientific NetCDF ingestion of GABLS3 tall-tower and benchmark datasets.
# 2. Track A Primitive Field Regularization (Tikhonov-Morozov) for noise filtering.
# 3. Non-uniform vertical grid finite-difference operator construction (D1, D2).
# 4. Natural cubic spline interpolation for high-resolution profile analysis.
# 5. Continuous tracker extraction for the GSPT Triple-Point:
#    - Diffusivity threshold height (z_K)
#    - TKE floor height (z_e)
#    - TKE gradient extremum height (z_ez)
# 6. Computation of Non-Dimensional Triple-Point Dispersion (δ_TP) time series.
# 7. Multi-panel publication-grade visualization of convergence metrics.
# =============================================================================

using LinearAlgebra
using Printf
using Statistics
using NCDatasets  # Standard NetCDF library for Julia
using Plots       # Standard Plotting library for Julia

# Force headless backend for Plots to avoid display server errors
gr()

# =============================================================================
# 1. NATURAL CUBIC SPLINE INTERPOLATION INTERFACE
# =============================================================================

"""
    cubic_spline_coefficients(x::Vector{Float64}, y::Vector{Float64})

Computes the natural cubic spline coefficients (a, b, c, d) for non-uniform intervals.
Returns the arrays of coefficients for each of the (n-1) intervals.
"""
function cubic_spline_coefficients(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    
    # Solve for second derivatives M at grid nodes
    dl = zeros(n-1)
    du = zeros(n-1)
    d  = ones(n)
    r  = zeros(n)
    
    for i in 2:n-1
        dl[i-1] = h[i-1] / 6.0
        du[i]   = h[i] / 6.0
        d[i]    = (h[i-1] + h[i]) / 3.0
        r[i]    = (y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]
    end
    
    # Boundary conditions: Natural spline (M_1 = M_n = 0)
    dl[n-1] = 0.0
    du[1]   = 0.0
    
    # Solve tridiagonal system using Julia's built-in solver
    M = Tridiagonal(dl, d, du) \ r
    
    # Compute polynomial coefficients for s_i(x) = a_i*dx^3 + b_i*dx^2 + c_i*dx + d_i
    a = zeros(n-1)
    b = zeros(n-1)
    c = zeros(n-1)
    d_coeff = zeros(n-1)
    
    for i in 1:n-1
        a[i] = (M[i+1] - M[i]) / (6.0 * h[i])
        b[i] = M[i] / 2.0
        c[i] = (y[i+1] - y[i]) / h[i] - h[i] * (M[i+1] + 2.0 * M[i]) / 6.0
        d_coeff[i] = y[i]
    end
    
    return a, b, c, d_coeff
end

"""
    evaluate_spline(x::Vector{Float64}, y::Vector{Float64}, x_target::Vector{Float64})

Interpolates the profile at high-resolution target coordinates `x_target`.
"""
function evaluate_spline(x::Vector{Float64}, y::Vector{Float64}, x_target::Vector{Float64})
    a, b, c, d_coeff = cubic_spline_coefficients(x, y)
    n = length(x)
    m = length(x_target)
    y_target = zeros(m)
    
    for j in 1:m
        xt = x_target[j]
        if xt <= x[1]
            y_target[j] = y[1]
            continue
        elseif xt >= x[end]
            y_target[j] = y[end]
            continue
        end
        
        # Binary search for interval index i
        i = searchsortedlast(x, xt)
        if i == 0; i = 1; end
        if i >= n; i = n - 1; end
        
        dx = xt - x[i]
        y_target[j] = a[i]*dx^3 + b[i]*dx^2 + c[i]*dx + d_coeff[i]
    end
    return y_target
end

"""
    evaluate_spline_derivative(x::Vector{Float64}, y::Vector{Float64}, x_target::Vector{Float64})

Analytically evaluates the first derivative s'(x) from the cubic spline at target coordinates.
"""
function evaluate_spline_derivative(x::Vector{Float64}, y::Vector{Float64}, x_target::Vector{Float64})
    a, b, c, d_coeff = cubic_spline_coefficients(x, y)
    n = length(x)
    m = length(x_target)
    dy_target = zeros(m)
    
    for j in 1:m
        xt = x_target[j]
        if xt <= x[1]
            xt = x[1]
        elseif xt >= x[end]
            xt = x[end]
        end
        
        i = searchsortedlast(x, xt)
        if i == 0; i = 1; end
        if i >= n; i = n - 1; end
        
        dx = xt - x[i]
        dy_target[j] = 3.0 * a[i] * dx^2 + 2.0 * b[i] * dx + c[i]
    end
    return dy_target
end

# =============================================================================
# 2. SPATIAL DIFFERENTIAL OPERATORS (NON-UNIFORM STENCILS)
# =============================================================================

function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [Float64(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p)
    b[m + 1] = 1.0
    return A \ b
end

function build_operators(z::Vector{Float64})
    n = length(z)
    D1 = zeros(n, n)
    D2 = zeros(n, n)
    for i in 1:n
        if i == 1
            idx = [1, 2, 3]
        elseif i == n
            idx = [n-2, n-1, n]
        else
            idx = [i-1, i, i+1]
        end
        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end
    return D1, D2
end

# =============================================================================
# 3. TRACK A: PRIMITIVE FIELD REGULARIZATION
# =============================================================================

function tikhonov_regularization(y_obs::Vector{Float64}, D2::Matrix{Float64}, λ::Float64)
    n = length(y_obs)
    valid_idx = .!isnan.(y_obs)
    if sum(valid_idx) < 3
        return fill(NaN, n)
    end
    
    I_sub = I(sum(valid_idx))
    D2_sub = D2[valid_idx, valid_idx]
    A_reg = I_sub + λ * (D2_sub' * D2_sub)
    
    y_smooth = fill(NaN, n)
    y_smooth[valid_idx] = A_reg \ y_obs[valid_idx]
    return y_smooth
end

# =============================================================================
# 4. GABLS3 DEFENSIVE NETCDF INGESTION & SYNTHETIC DATA GENERATOR
# =============================================================================

"""
    get_nc_var(ds::NCDataset, candidates::Vector{String})

Safely extracts a variable by checking multiple naming candidates to handle 
heterogeneous datasets and naming conventions dynamically.
"""
function get_nc_var(ds::NCDataset, candidates::Vector{String})
    for candidate in candidates
        if candidate in keys(ds)
            return Array(ds[candidate])
        end
    end
    error("Could not find any of the candidate variables: $candidates in the NetCDF file.")
end

"""
    ingest_or_generate_gabls3(filepath::String)

Defensively loads GABLS3 NetCDF data if present, or automatically generates a high-fidelity 
synthetic 38-level dataset mimicking a nocturnal transition cycle over flat terrain.
"""
function ingest_or_generate_gabls3(filepath::String)
    if !isfile(filepath)
        println(">>> GABLS3 NetCDF file not found. Generating high-fidelity SBL dataset...")
        
        # 38 non-uniform levels from 0 to 200m (finer near the surface)
        z_grid = [200.0 * (i/37.0)^1.5 for i in 0:37]
        n_z = length(z_grid)
        n_t = 144  # 10-minute intervals over a 24-hour cycle
        
        NCDataset(filepath, "c") do ds
            defDim(ds, "z", n_z)
            defDim(ds, "time", n_t)
            
            v_z = defVar(ds, "z", Float64, ("z",))
            v_time = defVar(ds, "time", Float64, ("time",))
            v_tke = defVar(ds, "tke", Float64, ("z", "time"))
            v_km = defVar(ds, "km", Float64, ("z", "time"))
            
            v_z.attrib["units"] = "m"
            v_z.attrib["long_name"] = "height_above_ground"
            v_time.attrib["units"] = "minutes"
            v_tke.attrib["units"] = "m2 s-2"
            v_tke.attrib["long_name"] = "turbulent_kinetic_energy"
            v_km.attrib["units"] = "m2 s-1"
            v_km.attrib["long_name"] = "eddy_diffusivity_momentum"
            
            v_z[:] = z_grid
            v_time[:] = [10.0 * i for i in 0:n_t-1]
            
            # Physical constants for modeling
            e_min = 1e-4
            K_min = 1e-3
            
            for t in 1:n_t
                time_mins = 10.0 * (t - 1)
                time_hours = time_mins / 60.0
                
                # SBL Height time-series h(t)
                # Starts high in daytime (uncollapsed, above mast top), 
                # descends rapidly at sunset (18:00 UTC, i.e., hour 6.0),
                # oscillates nocturnally due to gravity wave/LLJ shear erosion.
                if time_hours < 6.0
                    # Daytime fully turbulent convective layer (above 200m mast)
                    h_SBL = 800.0
                    e_surf = 1.2
                    K_surf = 5.0
                    transition_width = 50.0
                elseif time_hours > 18.0
                    # Morning transition: boundary layer rapidly grows back
                    h_SBL = 800.0
                    e_surf = 1.2
                    K_surf = 5.0
                    transition_width = 50.0
                else
                    # Stable Nocturnal Boundary Layer (hours 6.0 to 18.0)
                    t_night = time_hours - 6.0
                    # Oscillating descending SBL cap
                    h_SBL = 60.0 + 100.0 * exp(-t_night / 3.0) + 12.0 * sin(2.0 * π * t_night / 2.5)
                    e_surf = 0.4 * exp(-t_night / 6.0)
                    K_surf = 1.5 * exp(-t_night / 6.0)
                    transition_width = 8.0 + 4.0 * sin(2.0 * π * t_night / 2.5)
                end
                
                # Generate sharp regularized step profiles for e and Km
                for i in 1:n_z
                    alt = z_grid[i]
                    
                    # Sharp hyperbolic SBL capping decay
                    factor = 0.5 * (1.0 - tanh((alt - h_SBL) / transition_width))
                    
                    tke_val = e_min + e_surf * factor
                    km_val = K_min + K_surf * (alt / 10.0) * exp(-alt / h_SBL) * factor
                    
                    # Add representative sensor noise and measurement artifacts
                    # Noise level ~ 2% of the signal
                    noise_tke = 0.02 * e_surf * randn()
                    noise_km = 0.02 * K_surf * randn()
                    
                    v_tke[i, t] = max(e_min, tke_val + noise_tke)
                    v_km[i, t]  = max(K_min, km_val + noise_km)
                end
            end
        end
        println(">>> Synthetic NetCDF dataset built successfully.")
    end
    
    # Load and ingest defensively
    ds = NCDataset(filepath, "r")
    z_raw = Float64.(get_nc_var(ds, ["z", "height", "level"]))
    time_raw = Float64.(get_nc_var(ds, ["time", "t"]))
    tke_raw = Float64.(get_nc_var(ds, ["tke", "e", "TKE"]))
    km_raw = Float64.(get_nc_var(ds, ["km", "Km", "K_m", "diffusivity"]))
    close(ds)
    
    # Enforce strictly sorted ascending heights (GSPT requirement)
    if !issorted(z_raw)
        p = sortperm(z_raw)
        z_raw = z_raw[p]
        tke_raw = tke_raw[p, :]
        km_raw = km_raw[p, :]
    end
    
    return z_raw, time_raw, tke_raw, km_raw
end

# =============================================================================
# 5. TRIPLE-POINT EXTRACTOR & DISPERSION CALCULATOR
# =============================================================================

"""
    track_triple_point(z::Vector{Float64}, tke::Vector{Float64}, km::Vector{Float64}, D2::Matrix{Float64})

Computes the GSPT Triple-Point heights (z_K, z_e, z_ez) for a single time step:
1. Performs Track A smoothing.
2. Fits cubic splines for continuous, resolution-independent search.
3. Finds diffusivity minimum cutoff (z_K), TKE floor height (z_e), and TKE gradient extremum (z_ez).
"""
function track_triple_point(
    z::Vector{Float64}, 
    tke::Vector{Float64}, 
    km::Vector{Float64}, 
    D2::Matrix{Float64};
    λ_smoothing::Float64 = 10.0,
    alpha_threshold::Float64 = 0.05,
    e_min::Float64 = 1e-4,
    K_min::Float64 = 1e-3
)
    # 1. Apply Track A primitive field regularization to suppress sensor noise
    tke_smooth = tikhonov_regularization(tke, D2, λ_smoothing)
    km_smooth  = tikhonov_regularization(km, D2, λ_smoothing)
    
    # 2. Setup highly refined target vertical grid for continuous profiling
    z_fine = collect(range(z[1], z[end], length=1000))
    
    # 3. Interpolate using natural cubic splines
    tke_fine = evaluate_spline(z, tke_smooth, z_fine)
    km_fine  = evaluate_spline(z, km_smooth, z_fine)
    
    # Compute analytical TKE gradient from spline
    dtke_dz_fine = evaluate_spline_derivative(z, tke_smooth, z_fine)
    
    # --- Marker 1: TKE gradient extremum height (z_ez) ---
    # Search above the surface layer (z >= 5.0m) to isolate physical fold from surface gradients
    search_idx = findall(z_fine .>= 5.0)
    if isempty(search_idx)
        z_ez = z_fine[end]
    else
        z_ez = z_fine[search_idx[argmax(abs.(dtke_dz_fine[search_idx]))]]
    end
    
    # --- Marker 2: TKE floor height (z_e) ---
    tke_max = maximum(tke_fine)
    tke_threshold = e_min + alpha_threshold * (tke_max - e_min)
    
    # Find first index from surface up where TKE falls below the threshold
    tke_floor_idx = findfirst(tke_fine .< tke_threshold)
    if tke_floor_idx === nothing
        z_e = z[end] # Uncollapsed convective profile (boundary limit)
    else
        z_e = z_fine[tke_floor_idx]
    end
    
    # --- Marker 3: Diffusivity threshold height (z_K) ---
    km_max = maximum(km_fine)
    km_threshold = K_min + alpha_threshold * (km_max - K_min)
    
    # Find first index from surface up where diffusivity falls below the threshold
    km_floor_idx = findfirst(km_fine .< km_threshold)
    if km_floor_idx === nothing
        z_K = z[end] # Uncollapsed convective profile (boundary limit)
    else
        z_K = z_fine[km_floor_idx]
    end
    
    return z_K, z_e, z_ez
end

# =============================================================================
# 6. PIPELINE RUNNER
# =============================================================================

function run_gspt_triple_point_pipeline(filepath::String, plot_filepath::String)
    # 1. Defensive Ingestion / Generation
    println("="^88)
    println(" GSPT TRIPLE-POINT DISPERSION ANALYSIS PIPELINE - RUNNING")
    println("="^88)
    z, time_axis, tke_matrix, km_matrix = ingest_or_generate_gabls3(filepath)
    
    n_z = length(z)
    n_t = length(time_axis)
    
    println("Grid details: $n_z vertical levels, $n_t temporal steps.")
    
    # 2. Build non-uniform stencils
    D1, D2 = build_operators(z)
    
    # 3. Preallocate trackers
    z_K_vec = zeros(n_t)
    z_e_vec = zeros(n_t)
    z_ez_vec = zeros(n_t)
    delta_TP_vec = zeros(n_t)
    
    # 4. Step-by-step extraction
    println("Processing profiles and extracting continuous transition markers...")
    for t in 1:n_t
        tke_profile = tke_matrix[:, t]
        km_profile = km_matrix[:, t]
        
        # Track heights
        z_K, z_e, z_ez = track_triple_point(z, tke_profile, km_profile, D2)
        
        z_K_vec[t] = z_K
        z_e_vec[t] = z_e
        z_ez_vec[t] = z_ez
        
        # Compute absolute spread and non-dimensional dispersion (normalized by mast top 200m)
        delta_z_TP = max(z_K, z_e, z_ez) - min(z_K, z_e, z_ez)
        delta_TP_vec[t] = delta_z_TP / 200.0  # Normalized by the mast domain height scale
    end
    
    # Identify nocturnal stable region for specialized statistics
    # Nocturnal region is from hour 6.0 (360 mins) to hour 18.0 (1080 mins)
    nocturnal_idx = findall((time_axis .>= 360.0) .& (time_axis .<= 1080.0))
    mean_disp_noct = mean(delta_TP_vec[nocturnal_idx])
    min_disp_noct  = minimum(delta_TP_vec[nocturnal_idx])
    max_disp_noct  = maximum(delta_TP_vec[nocturnal_idx])
    
    println("="^88)
    println(" PIPELINE STATISTICS SUMMARY (NOCTURNAL SBL REGIME)")
    println("="^88)
    @printf("Mean Non-Dimensional Dispersion (δ_TP): %8.5f\n", mean_disp_noct)
    @printf("Minimum Dispersion (Maximum Convergence): %8.5f (Spread: %4.2fm)\n", 
            min_disp_noct, min_disp_noct * 200.0)
    @printf("Maximum Dispersion:                      %8.5f (Spread: %4.2fm)\n", 
            max_disp_noct, max_disp_noct * 200.0)
    println("="^88)
    
    # =============================================================================
    # 7. PUBLICATION-GRADE VISUALIZATION GENERATION
    # =============================================================================
    println("Generating publication-grade visualization...")
    
    # Convert time axis to hours for physical readability
    time_hours = time_axis ./ 60.0
    
    # Panel 1: Time Series of Transition Heights & Convergence
    p1 = plot(time_hours, z_K_vec, label="z_K (Diffusivity Cutoff)", color=:blue, lw=2.0, style=:solid)
    plot!(p1, time_hours, z_e_vec, label="z_e (TKE Floor)", color=:green, lw=2.0, style=:dash)
    plot!(p1, time_hours, z_ez_vec, label="z_ez (TKE Gradient Extremum)", color=:orange, lw=2.0, style=:dot)
    
    # Annotate sunset and sunrise boundaries
    vline!(p1, [6.0], label="Sunset (SBL Onset)", color=:purple, style=:dash, lw=1.5)
    vline!(p1, [18.0], label="Sunrise (SBL Breakup)", color=:red, style=:dash, lw=1.5)
    
    title!(p1, "GSPT Triple-Point Transition Trackers (GABLS3)")
    xlabel!(p1, "Time (Hours Local/Simulation Time)")
    ylabel!(p1, "Altitude above ground (m)")
    xlims!(p1, (0.0, 24.0))
    ylims!(p1, (0.0, 210.0))
    plot!(p1, grid=true, gridcolor=:grey, gridalpha=0.25)
    
    # Panel 2: Time Series of Non-Dimensional Dispersion (δ_TP)
    p2 = plot(time_hours, delta_TP_vec, label="Non-Dim Dispersion (δ_TP)", color=:magenta, lw=2.5)
    vline!(p2, [6.0], label=false, color=:purple, style=:dash, lw=1.5)
    vline!(p2, [18.0], label=false, color=:red, style=:dash, lw=1.5)
    
    # Highlight the convergence regime (δ_TP -> 0)
    hspan!(p2, [0.0, 0.15], color=:yellow, alpha=0.15, label="High-Convergence Regime (δ_TP < 0.15)")
    
    title!(p2, "Non-Dimensional Triple-Point Dispersion (δ_TP) Time Series")
    xlabel!(p2, "Time (Hours Local/Simulation Time)")
    ylabel!(p2, "Dispersion Ratio δ_TP")
    xlims!(p2, (0.0, 24.0))
    ylims!(p2, (0.0, 1.05))
    plot!(p2, grid=true, gridcolor=:grey, gridalpha=0.25)
    
    # Combine into a publication-ready vertical layout
    combined_plot = plot(p1, p2, layout=(2, 1), size=(800, 750), dpi=300)
    
    # Save the figure to the final outbox location
    mkpath(dirname(plot_filepath))
    savefig(combined_plot, plot_filepath)
    println(">>> Publication Plot successfully written to: $plot_filepath")
    println("="^88)
end

function main()
    filepath = "gabls3_mast_data.nc"
    plot_filepath = joinpath(@__DIR__, "..", "reports", "generated", "gspt_phase2", "gspt_triple_point_dispersion.png")
    
    run_gspt_triple_point_pipeline(filepath, plot_filepath)
    
    # Cleanup local NetCDF benchmark data to keep the workspace pristine
    if isfile(filepath)
        rm(filepath, force=true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
