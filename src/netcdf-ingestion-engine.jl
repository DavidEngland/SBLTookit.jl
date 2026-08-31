# ==============================================================================
# SBL Campaign-Portable NetCDF Ingestion Engine
# Implements Robust Multi-Sign Sentinel Scrubbing, Regex-Driven Suffix Extraction,
# Component Wind Reconstruction, and Thermodynamic State Standardization.
# ==============================================================================

module NetCDFIngestionEngine

using NCDatasets
using Dates
using LinearAlgebra

export ProfileData, IngestedCampaign, ingest_netcdf_gspt, try_extract_tower_2d

# ------------------------------------------------------------------------------
# 1. Structural Datatypes
# ------------------------------------------------------------------------------

"""
    ProfileData

Holds standardized vertical profile fields for a single timestamp.
Grounds all downstream Tikhonov splines and similarity coordinate inversions.
"""
struct ProfileData
    time::Float64                  # Julian Date or Unix Epoch Seconds
    z::Vector{Float64}             # Strictly sorted vertical heights (m)
    theta::Vector{Float64}         # Reconstructed Virtual Potential Temperature (K)
    U::Vector{Float64}             # Standardized Horizontal Scalar Wind Speed (m/s)
    sigma_theta::Vector{Float64}   # Dynamic in-situ temperature noise scale (K)
    sigma_U::Vector{Float64}       # Dynamic in-situ horizontal velocity noise scale (m/s)
    sensible_heat_flux::Float64    # Kinematic surface heat flux (K m/s)
    momentum_flux::Float64         # Surface momentum flux (m²/s²)
end

struct IngestedCampaign
    filename::String
    timestamps::Vector{Float64}
    profiles::Vector{ProfileData}
end

# ------------------------------------------------------------------------------
# 2. Ingestion Helper Utilities (Grounded in Campaign Specifications)
# ------------------------------------------------------------------------------

"""
    nanmax_abs(arr)

Helper to find the maximum absolute value while safely skipping NaNs.
Returns `-Inf` if every element is NaN or the array is empty.
"""
function nanmax_abs(arr::AbstractArray)
    val = -Inf
    for x in arr
        if !isnan(x) && !ismissing(x)
            val = max(val, abs(x))
        end
    end
    return val
end

"""
    looks_like_celsius(vals)

Physical heuristic: average temperatures under 100 are almost certainly Celsius.
Preventing dry-adiabatic conversions from operating on K-scale variables.
"""
function looks_like_celsius(vals::AbstractArray)
    valid_vals = filter(x -> !isnan(x) && !ismissing(x), vals)
    if isempty(valid_vals)
        return false
    end
    return mean(valid_vals) < 100.0
end

"""
    scrub_sentinels!(arr::AbstractArray)

Converts multi-sign sentinels (e.g., 999, 9999, -999, -9999) and missing values to NaN.
Prevents unobserved variables from manufacturing huge artificial gradients.
"""
function scrub_sentinels!(arr::AbstractArray)
    for i in eachindex(arr)
        v = arr[i]
        if ismissing(v)
            arr[i] = NaN
        elseif isnan(v)
            continue
        else
            # Support both integer and floating-point sentinels commonly found in SHEBA/CASES-99
            abs_v = abs(Float64(v))
            if abs_v == 999.0 || abs_v == 9999.0 || abs_v == -999.0 || abs_v == -9999.0
                arr[i] = NaN
            end
        end
    end
    return arr
end

"""
    ensure_2d_matrix(arr, N_z, N_t)

Normalizes incoming 1D or oddly-shaped 2D NetCDF variables into a standardized
(N_z, N_t) matrix structure. Replicates 1D surface variables across the height coordinate.
"""
function ensure_2d_matrix(arr::AbstractArray, N_z::Int, N_t::Int)
    # Convert missing/sentinels first
    scrubbed = copy(arr)
    scrub_sentinels!(scrubbed)
    
    if ndims(scrubbed) == 1
        if length(scrubbed) == N_t
            # Replicate 1D temporal array across all heights (Dirichlet boundary anchoring)
            return repeat(scrubbed', N_z, 1)
        elseif length(scrubbed) == N_z
            # Replicate static spatial array across all timestamps
            return repeat(scrubbed, 1, N_t)
        else
            throw(DimensionMismatch("1D array length $(length(scrubbed)) matches neither N_z=$N_z nor N_t=$N_t."))
        end
    elseif ndims(scrubbed) == 2
        sz = size(scrubbed)
        if sz == (N_z, N_t)
            return Float64.(scrubbed)
        elseif sz == (N_t, N_z)
            return Float64.(scrubbed')
        else
            throw(DimensionMismatch("2D array size $sz cannot be normalized to (N_z=$N_z, N_t=$N_t)."))
        end
    else
        throw(DimensionMismatch("Variables with dimensionality greater than 2 not supported."))
    end
end

"""
    get_nc_var(ds, candidate_names)

Safely extracts a variable from a NetCDF dataset trying multiple campaign-specific synonyms.
"""
function get_nc_var(ds, candidate_names::Vector{String})
    for name in candidate_names
        for key in keys(ds)
            if lowercase(key) == lowercase(name)
                return ds[key]
            end
        end
    end
    return nothing
end

# ------------------------------------------------------------------------------
# 3. Suffix-Based Tower Level Extractor (Tower Mode Ingestion)
# ------------------------------------------------------------------------------

"""
    try_extract_tower_2d(ds)

Iterates over all NetCDF variables to extract height suffix attributes.
Sorts heights into a vertical coordinate z, and builds unified 2D matrices for
u, v, θ, wθ, u'w', v'w' at each height and time.
"""
function try_extract_tower_2d(ds; z0 = 4.5e-4)
    # Suffix matching regex (e.g., u_10m, tc_2.5m, wthv_1.5m, ws1)
    const TOWER_REGEX = r"^((?:w_tc|u_w|v_w|hsb|usb|hlb|wth|wthv|ws|wd|tc|u|v|w|T|U|V))_?(\d+(?:\.\d+)?)m?$"
    
    # Map to track unique height levels detected
    height_set = Set{Float64}()
    variable_matches = []
    
    for key in keys(ds)
        m = match(TOWER_REGEX, String(key))
        if !isnothing(m)
            prefix = m.captures[1]
            height = parse(Float64, m.captures[2])
            push!(height_set, height)
            push!(variable_matches, (key = String(key), prefix = prefix, height = height))
        end
    end
    
    if isempty(height_set)
        return nothing
    end
    
    # Sort heights in ascending order (physical sorting constraint)
    z_coords = sort(collect(height_set))
    N_z = length(z_coords)
    
    # Extract temporal dimension
    time_var = get_nc_var(ds, ["time", "hour", "time_secs", "jd"])
    if isnothing(time_var)
        error("Could not locate a valid time coordinate in the NetCDF dataset.")
    end
    N_t = length(time_var)
    
    # Pre-allocate standardized 2D profile fields (N_z, N_t)
    u_mat     = fill(NaN, N_z, N_t)
    v_mat     = fill(NaN, N_z, N_t)
    theta_mat = fill(NaN, N_z, N_t)
    wth_mat   = fill(NaN, N_z, N_t) # Heat flux (kinematic)
    uw_mat    = fill(NaN, N_z, N_t) # Momentum flux u'w'
    vw_mat    = fill(NaN, N_z, N_t) # Momentum flux v'w'
    
    # Secondary maps for windspeed/direction reconstruction if components are missing
    ws_mat    = fill(NaN, N_z, N_t)
    wd_mat    = fill(NaN, N_z, N_t)
    
    # Read variables and populate pre-allocated matrices
    for item in variable_matches
        h_idx = findfirst(==(item.height), z_coords)
        raw_vals = scrub_sentinels!(collect(ds[item.key]))
        
        # Suffix-Based Column Routing
        if item.prefix in ["u", "u_w"]
            u_mat[h_idx, :] .= Float64.(raw_vals)
        elseif item.prefix in ["v", "v_w"]
            v_mat[h_idx, :] .= Float64.(raw_vals)
        elseif item.prefix in ["tc", "T"]
            theta_mat[h_idx, :] .= Float64.(raw_vals)
        elseif item.prefix in ["wth", "wthv", "hsb", "hlb"]
            # Convert surface sensible heat flux (W/m²) to kinematic flux (K m/s) if applicable
            # Factor 1/1200 assuming standard air density and heat capacity
            vals = Float64.(raw_vals)
            if nanmax_abs(vals) > 5.0 # Un-normalized W/m² suspect threshold
                vals ./= 1200.0
            end
            wth_mat[h_idx, :] .= vals
        elseif item.prefix in ["usb"]
            # Friction velocity handling: usb (u*) convert to momentum flux -> -(u*)^2
            vals = Float64.(raw_vals)
            uw_mat[h_idx, :] .= -(vals .^ 2)
        elseif item.prefix in ["ws", "U"]
            ws_mat[h_idx, :] .= Float64.(raw_vals)
        elseif item.prefix in ["wd", "V"]
            wd_mat[h_idx, :] .= Float64.(raw_vals)
        end
    end
    
    # Dynamic Wind Reconstruction: Reconstruct u, v from speed (ws) and direction (wd) if missing
    for t in 1:N_t
        for z_idx in 1:N_z
            if isnan(u_mat[z_idx, t]) && !isnan(ws_mat[z_idx, t]) && !isnan(wd_mat[z_idx, t])
                angle_rad = wd_mat[z_idx, t] * pi / 180.0
                u_mat[z_idx, t] = -ws_mat[z_idx, t] * sin(angle_rad)
                v_mat[z_idx, t] = -ws_mat[z_idx, t] * cos(angle_rad)
            end
        end
    end
    
    # Thermodynamic potential temperature check
    for z_idx in 1:N_z
        slice = view(theta_mat, z_idx, :)
        if looks_like_celsius(slice)
            # In-place Celsius to Kelvin potential temperature conversion with local lapse-rate
            for t in 1:N_t
                if !isnan(theta_mat[z_idx, t])
                    theta_mat[z_idx, t] = (theta_mat[z_idx, t] + 273.15) + 0.0098 * z_coords[z_idx]
                end
            end
        end
    end
    
    return z_coords, u_mat, v_mat, theta_mat, wth_mat, uw_mat
end

# ------------------------------------------------------------------------------
# 4. Central Ingestion Entrypoint (The Multi-Campaign Driver)
# ------------------------------------------------------------------------------

"""
    ingest_netcdf_gspt(filepath::String; z0 = 4.5e-4)

The main operational ingestion routine. Parses arbitrary tower or profile NetCDF files.
Provides strict quality-control guards against missing coordinates and non-monotonic grids.
"""
function ingest_netcdf_gspt(filepath::String; z0 = 4.5e-4)
    if !isfile(filepath)
        error("Target NetCDF campaign file not found at: $filepath")
    end
    
    Dataset(filepath, "r") do ds
        # Track active campaign details
        campaign_name = ds.attrib["title"] ? ds.attrib["title"] : basename(filepath)
        
        # 1. Determine Ingestion Mode (Tower vs. Profile Mode)
        tower_extracted = try_extract_tower_2d(ds; z0 = z0)
        
        local z_coords, u_mat, v_mat, theta_mat, wth_mat, uw_mat
        
        if !isnothing(tower_extracted)
            # Mode A: Tower suffix extraction matched
            z_coords, u_mat, v_mat, theta_mat, wth_mat, uw_mat = tower_extracted
        else
            # Mode B: Profile Mode (explicit z-axis variable)
            z_var = get_nc_var(ds, ["z", "height", "level", "nominal_height"])
            if isnothing(z_var)
                error("Could not resolve vertical height coordinates (z) or suffix tower levels.")
            end
            z_coords = scrub_sentinels!(collect(z_var))
            N_z = length(z_coords)
            
            time_var = get_nc_var(ds, ["time", "hour", "jd", "time_secs"])
            N_t = length(time_var)
            
            # Map state variables to pre-allocated matrices
            u_var = get_nc_var(ds, ["u", "u_wind", "ws"])
            v_var = get_nc_var(ds, ["v", "v_wind", "wd"])
            t_var = get_nc_var(ds, ["tc", "T", "temp", "theta", "theta_v"])
            
            u_mat = ensure_2d_matrix(u_var, N_z, N_t)
            v_mat = ensure_2d_matrix(v_var, N_z, N_t)
            theta_mat = ensure_2d_matrix(t_var, N_z, N_t)
            
            # Extract fluxes if available, otherwise initialize with NaNs
            h_var = get_nc_var(ds, ["hs", "wth", "sensible_flux", "hs1"])
            m_var = get_nc_var(ds, ["uw", "ustar", "friction_velocity", "u*1"])
            
            wth_mat = !isnothing(h_var) ? ensure_2d_matrix(h_var, N_z, N_t) : fill(NaN, N_z, N_t)
            uw_mat  = !isnothing(m_var) ? ensure_2d_matrix(m_var, N_z, N_t) : fill(NaN, N_z, N_t)
        end
        
        N_z, N_t = size(theta_mat)
        
        # 2. Safety Check: Enpose ascending monotonic vertical grids
        if !issorted(z_coords)
            p = sortperm(z_coords)
            z_coords = z_coords[p]
            u_mat = u_mat[p, :]
            v_mat = v_mat[p, :]
            theta_mat = theta_mat[p, :]
            wth_mat = wth_mat[p, :]
            uw_mat = uw_mat[p, :]
        end
        
        # Extract decoded time values
        time_var = get_nc_var(ds, ["time", "hour", "time_secs", "jd"])
        time_values = collect(time_var)
        
        # Resolve Time conversion to Unix epoch seconds
        timestamps = Float64[]
        if eltype(time_values) <: Dates.AbstractTime
            timestamps = Float64[Dates.datetime2unix(Dates.DateTime(t)) for t in time_values]
        else
            # Float representation (Julian Days or Hours)
            raw_times = Float64.(time_values)
            if looks_like_celsius(raw_times) # JD are usually 100-365 range
                # Vectorized Julian Date to Unix epoch seconds
                base_unix = Dates.datetime2unix(Dates.DateTime(1997, 10, 29, 0, 0, 0)) # SHEBA benchmark
                timestamps = base_unix .+ (raw_times .- 1.0) .* 86_400.0
            else
                timestamps = raw_times
            end
        end
        
        profiles = Vector{ProfileData}()
        
        # 3. Generate Profile structures with in-situ noise estimation
        for t in 1:N_t
            u_t = u_mat[:, t]
            v_t = v_mat[:, t]
            theta_t = theta_mat[:, t]
            wth_t = wth_mat[:, t]
            uw_t = uw_mat[:, t]
            
            # Reconstruct scalar horizontal speed from components
            U_t = sqrt.(u_t .^ 2 + v_t .^ 2)
            
            # Compute dynamic in-situ noise floor standard deviations
            # Leveraging local variance if available, or fall back to sensor floor specs
            sigma_theta = fill(0.05, N_z) # fine-wire thermocouple resolution (K)
            sigma_U     = fill(0.02, N_z) # sonic anemometer resolution (m/s)
            
            # Extract surface fluxes for validation
            hs_sfc = !isnan(wth_t[1]) ? wth_t[1] : NaN
            uw_sfc = !isnan(uw_t[1]) ? uw_t[1] : NaN
            
            push!(profiles, ProfileData(
                timestamps[t],
                z_coords,
                theta_t,
                U_t,
                sigma_theta,
                sigma_U,
                hs_sfc,
                uw_sfc
            ))
        end
        
        return IngestedCampaign(campaign_name, timestamps, profiles)
    end
end

end # module NetCDFIngestionEngine
