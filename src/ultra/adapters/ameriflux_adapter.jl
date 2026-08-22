module AmeriFluxAdapter

using Dates
using DataFrames
using CSV
using JSON3
using ..CoreTypes

export load_ameriflux_registry, get_site_metadata,
    extract_ameriflux_observations, extract_ameriflux_multi_site

# Path to default station registry
const AMERIFLUX_STATIONS_PATH = joinpath(@__DIR__, "..", "..", "..", "data", "ameriflux", "stations.json")

"""
    load_ameriflux_registry(registry_path=AMERIFLUX_STATIONS_PATH) -> Dict

Load AmeriFlux station metadata from JSON registry.
"""
function load_ameriflux_registry(registry_path::String=AMERIFLUX_STATIONS_PATH)::Dict
    if !isfile(registry_path)
        @warn "AmeriFlux registry not found at $registry_path; using fallback defaults"
        return Dict()
    end

    try
        registry_json = open(registry_path) do f
            JSON3.read(f, Dict)
        end
        registry = Dict()

        for site in get(registry_json, "sites", [])
            site_id = site["id"]
            registry[site_id] = (
                heights=Float64.(get(site, "measurement_heights_m", Float64[])),
                z0m=Float64(get(site, "z0m", 0.1)),
                d_displacement=Float64(get(site, "d_displacement", 0.0))
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

"""
    get_site_metadata(site_id::String)

Retrieve site-specific heights and roughness from registry.
"""
function get_site_metadata(site_id::String)
    if haskey(AMERIFLUX_REGISTRY, site_id)
        return AMERIFLUX_REGISTRY[site_id]
    else
        @warn "Site $site_id not in registry; using generic defaults"
        return (
            heights=[2.0, 10.0, 20.0, 40.0],
            z0m=0.1,
            d_displacement=0.0
        )
    end
end

"""
    extract_ameriflux_observations(site_id::String, csv_path::String; min_levels::Int=3)

Extract StandardizedBLObservation records from an AmeriFlux BASE CSV file.
"""
function extract_ameriflux_observations(site_id::String, csv_path::String; min_levels::Int=3)
    site_meta = get_site_metadata(site_id)
    site_heights = site_meta.heights
    z0m = site_meta.z0m

    df = CSV.read(csv_path, DataFrame)
    observations = StandardizedBLObservation[]

    for row in eachrow(df)
        timestamp_val = row["TIMESTAMP_START"]
        if ismissing(timestamp_val) || timestamp_val == -9999
            continue
        end

        timestamp_str = string(Int(timestamp_val))

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
            continue
        end

        temp_values = Float64[]
        heights_m = Float64[]

        for (h_idx, h) in enumerate(site_heights)
            col_name_v1 = "TA_$(h_idx)_1_1"
            col_name_v2 = "TA_1_$(h_idx)_1"

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
                push!(temp_values, val)
            end
        end

        n_valid = length(temp_values)
        if n_valid < min_levels
            continue
        end

        ustar = 0.3
        if haskey(row, "USTAR") && !ismissing(row["USTAR"])
            ustar_val = row["USTAR"]
            if ustar_val != -9999 && !isnan(ustar_val)
                ustar = max(0.01, Float64(ustar_val))
            end
        end

        L_obukhov = 9999.0
        if haskey(row, "ZL") && !ismissing(row["ZL"])
            zl = row["ZL"]
            if zl != -9999 && !isnan(zl) && abs(zl) > 1e-5
                try
                    L_calc = heights_m[1] / Float64(zl)
                    L_obukhov = max(-9999.0, min(9999.0, L_calc))
                catch
                    L_obukhov = 9999.0
                end
            end
        end

        obs = StandardizedBLObservation(
            datetime=dt,
            campaign="AMERIFLUX",
            heights=heights_m,
            values=temp_values,
            ustar=ustar,
            L_obukhov=L_obukhov,
            z0m=z0m,
            robust_for_eta3=(n_valid >= 3),
            n_valid_levels=n_valid
        )

        push!(observations, obs)
    end

    return observations
end

"""
    extract_ameriflux_multi_site(site_dir="data/ameriflux/base_badm_extracted"; target_sites=nothing)

Extract observations from all AmeriFlux sites in a directory.
"""
function extract_ameriflux_multi_site(site_dir::String="data/ameriflux/base_badm_extracted";
    target_sites::Union{Nothing,Vector{String}}=nothing)
    rows = []

    if !isdir(site_dir)
        @warn "Directory $site_dir does not exist"
        return DataFrame()
    end

    site_dirs = readdir(site_dir)
    if !isnothing(target_sites)
        site_dirs = filter(s -> s in target_sites, site_dirs)
    end

    for site_id in site_dirs
        site_path = joinpath(site_dir, site_id)
        if !isdir(site_path)
            continue
        end

        csv_files = filter(f -> endswith(f, ".csv"), readdir(site_path))
        if isempty(csv_files)
            @warn "No CSV found for $site_id in $site_path"
            continue
        end

        csv_path = joinpath(site_path, csv_files[1])

        try
            obs_list = extract_ameriflux_observations(site_id, csv_path)
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

            @info "Loaded $(length(obs_list)) observations from $site_id"
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

end # module AmeriFluxAdapter