"""
AmeriFluxAdapter.jl
===================
Adapter for AmeriFlux BASE-BADM data products (half-hourly CSV).

Converts half-hourly tower measurements from AmeriFlux FP Standard format
into StandardizedBLObservation contracts for cross-campaign ingestion.

AmeriFlux sites: 60+ NEON + research tower network across North America.
Data: https://ameriflux.lbl.gov/ (open science, CC-BY-4.0 license)

FP Standard Variable Naming:
  - TA_H_V_R: Air temperature at height H, vertical position V, replication R
  - RH_H_V_R: Relative humidity
  - WS_H_V_R: Wind speed
  - WD_H_V_R: Wind direction
  - USTAR: Friction velocity (m/s)
  - ZL: Stability parameter z/L (dimensionless, Monin-Obukhov)
  - PA: Atmospheric pressure (kPa)

Missing data indicator: -9999 (set to NaN)

CRITICAL: Heights are site-specific and loaded from data/ameriflux/stations.json
to ensure correct spatial gridding of eta operators across heterogeneous tower configs.
"""

# Load station registry at module initialization
const AMERIFLUX_STATIONS_PATH = "data/ameriflux/stations.json"

function load_ameriflux_registry(registry_path::String=AMERIFLUX_STATIONS_PATH)::Dict
    """
    Load AmeriFlux station metadata from JSON registry.
    
    Returns: Dict mapping site_id → (heights, z0m, d_displacement)
    """
    if !isfile(registry_path)
        @warn "AmeriFlux registry not found at $registry_path; using fallback defaults"
        return Dict()
    end
    
    try
        # Read JSON file using JSON3 (standard in this project)
        registry_json = open(registry_path) do f
            JSON3.read(f, Dict)
        end
        registry = Dict()
        
        for site in get(registry_json, "sites", [])
            site_id = site["id"]
            registry[site_id] = (
                heights = get(site, "measurement_heights_m", Float64[]),
                z0m = get(site, "z0m", 0.1),
                d_displacement = get(site, "d_displacement", 0.0)
            )
        end
        
        return registry
    catch e
        @warn "Failed to parse AmeriFlux registry: $e; using fallback defaults"
        return Dict()
    end
end

# Lazy-load registry
const AMERIFLUX_REGISTRY = load_ameriflux_registry()

function get_site_metadata(site_id::String)
    """
    Retrieve site-specific heights and roughness from registry.
    
    Returns: (heights_m, z0m, d_displacement) or defaults if not found
    """
    if haskey(AMERIFLUX_REGISTRY, site_id)
        return AMERIFLUX_REGISTRY[site_id]
    else
        # Fallback for sites not in registry
        @warn "Site $site_id not in registry; using generic defaults"
        return (
            heights = [2.0, 10.0, 20.0, 40.0],
            z0m = 0.1,
            d_displacement = 0.0
        )
    end
end


function extract_ameriflux_observations(site_id::String, csv_path::String;
                                        min_levels::Int=3)
    """
    Extract StandardizedBLObservation records from AmeriFlux BASE CSV.
    
    CRITICAL: Uses site-specific measurement heights from stations.json registry
    to correctly map TA_1_1_1, TA_2_1_1, etc. to physical heights.
    
    Args:
        site_id: AmeriFlux site code (e.g., "US-ARM", "US-xCP")
        csv_path: Path to extracted CSV file
        min_levels: Minimum number of valid heights required for a profile (default: 3)
    
    Returns:
        Vector{StandardizedBLObservation}
    
    Raises:
        @warn if site not found in registry; continues with fallback heights
    """
    
    # Retrieve site-specific configuration from registry
    site_meta = get_site_metadata(site_id)
    site_heights = site_meta.heights
    z0m = site_meta.z0m
    d_displacement = site_meta.d_displacement
    
    # Read CSV
    df = CSV.read(csv_path, DataFrame)
    
    observations = StandardizedBLObservation[]
    
    for row in eachrow(df)
        # Extract timestamp (TIMESTAMP_START format: YYYYMMDDHHMM)
        timestamp_val = row["TIMESTAMP_START"]
        if ismissing(timestamp_val) || timestamp_val == -9999
            continue
        end
        
        timestamp_str = string(Int(timestamp_val))
        
        # Parse YYYYMMDDHHMM → DateTime
        dt = try
            if length(timestamp_str) >= 12
                date_part = timestamp_str[1:8]
                hour_part = parse(Int, timestamp_str[9:10])
                min_part = parse(Int, timestamp_str[11:12])
                DateTime(date_part, DateFormat("yyyymmdd")) + Hour(hour_part) + Minute(min_part)
            else
                continue
            end
        catch
            continue  # Skip malformed timestamps
        end
        
        # DYNAMIC HEIGHT PARSING: Use site-specific heights from registry
        # For each height in the registry, look for the corresponding column
        # The column index h_idx is sequential (1, 2, 3, ...) matching the registry order
        
        temp_values = Float64[]
        heights_m = Float64[]
        
        for (h_idx, h) in enumerate(site_heights)
            # Try both FP Standard notations (some towers use TA_H_V_R vs TA_1_H_R)
            col_name_v1 = "TA_$(h_idx)_1_1"  # Standard: index by height position
            col_name_v2 = "TA_1_$(h_idx)_1"  # Alternative: vertical tier notation
            
            val = nothing
            for col_name in [col_name_v1, col_name_v2]
                if haskey(row, col_name) && !ismissing(row[col_name])
                    raw_val = row[col_name]
                    if raw_val != -9999 && !isnan(raw_val)
                        val = Float64(raw_val)
                        break
                    end
                end
            end
            
            if !isnothing(val)
                push!(heights_m, h)
                push!(temp_values, val)  # Temperature in Celsius
            end
        end
        
        # Skip profiles with insufficient levels
        n_valid = length(temp_values)
        if n_valid < min_levels
            continue
        end
        
        # Extract stability and friction parameters with defensive clamping
        
        # USTAR: friction velocity (m/s)
        ustar = 0.3  # Default fallback (m/s)
        if haskey(row, "USTAR") && !ismissing(row["USTAR"])
            ustar_val = row["USTAR"]
            if ustar_val != -9999 && !isnan(ustar_val)
                ustar = max(0.01, Float64(ustar_val))  # Clamp to avoid zero
            end
        end
        
        # DEFENSIVE CLAMP: L_obukhov boundary
        # Under neutral conditions (z/L → 0), L_obukhov → Inf, which breaks LaTeX plotting.
        # Use 9999 m as physical limit (effectively neutral for all practical heights).
        L_obukhov = 9999.0  # Default: effectively neutral
        if haskey(row, "ZL") && !ismissing(row["ZL"])
            zl = row["ZL"]
            if zl != -9999 && !isnan(zl) && abs(zl) > 1e-5  # Avoid division by near-zero
                try
                    # L = z / (z/L) where z is first measurement height
                    L_calc = heights_m[1] / Float64(zl)
                    # Clamp to ±9999 m (practical boundary for neutral/near-neutral)
                    L_obukhov = max(-9999.0, min(9999.0, L_calc))
                catch
                    L_obukhov = 9999.0  # On error, use neutral default
                end
            end
        end
        
        # Create StandardizedBLObservation with site-specific parameters
        obs = StandardizedBLObservation(
            datetime=dt,
            campaign="AMERIFLUX",
            heights=heights_m,
            values=temp_values,
            ustar=ustar,
            L_obukhov=L_obukhov,
            z0m=z0m,  # Site-specific roughness from registry
            robust_for_eta3=(n_valid >= 3),
            n_valid_levels=n_valid
        )
        
        push!(observations, obs)
    end
    
    return observations
end


function extract_ameriflux_multi_site(site_dir::String="data/ameriflux/base_badm_extracted";
                                       target_sites::Union{Nothing, Vector{String}}=nothing)
    """
    Extract observations from all AmeriFlux sites in a directory.
    
    Dynamically walks site_dir, loads height configs from registry, and processes
    each site's CSV independently. Handles missing registry entries gracefully.
    
    Args:
        site_dir: Root directory containing subdirectories like US-ARM/, US-xCP/, etc.
        target_sites: If provided, only process these site IDs. Otherwise, process all.
    
    Returns:
        DataFrame with columns: site_id, datetime, campaign, n_valid_levels, robust_for_eta3, ustar, L_obukhov, z0m
    """
    
    rows = []
    
    site_dirs = readdir(site_dir)
    if !isnothing(target_sites)
        site_dirs = filter(s -> s in target_sites, site_dirs)
    end
    
    for site_id in site_dirs
        site_path = joinpath(site_dir, site_id)
        if !isdir(site_path)
            continue
        end
        
        # Find CSV file (wildcard match for flexibility in naming)
        csv_files = filter(f -> endswith(f, ".csv"), readdir(site_path))
        if isempty(csv_files)
            @warn "No CSV found for $site_id in $site_path"
            continue
        end
        
        csv_path = joinpath(site_path, csv_files[1])
        
        try
            obs_list = extract_ameriflux_observations(site_id, csv_path)
            
            # Get site-specific roughness for DataFrame
            site_meta = get_site_metadata(site_id)
            z0m_val = site_meta.z0m
            
            for obs in obs_list
                push!(rows, (
                    site_id=site_id,
                    datetime=obs.datetime,
                    campaign=obs.campaign,
                    n_valid_levels=obs.n_valid_levels,
                    robust_for_eta3=obs.robust_for_eta3,
                    ustar=obs.ustar,
                    L_obukhov=obs.L_obukhov,
                    z0m=z0m_val
                ))
            end
            
            @info "Loaded $(length(obs_list)) observations from $site_id ($(length(obs_list) > 0 ? "$(size(first(obs_list).heights)[1]) levels" : "N/A"))"
        catch e
            @warn "Failed to load $site_id: $(sprint(showerror, e))"
            continue
        end
    end
    
    if isempty(rows)
        @warn "No observations extracted from any sites in $site_dir"
        return DataFrame()
    end
    
    return DataFrame(rows)
end
