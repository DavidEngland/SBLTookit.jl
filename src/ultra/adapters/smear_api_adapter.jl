module SmearApiAdapter

using HTTP
using JSON3
using Dates
using DataFrames
using ..CoreTypes

export SmearStationConfig, fetch_smear_observations, parse_smear_json

struct SmearStationConfig
    station_id::Int
    name::String
    z0m::Float64
    reference_height::Float64
    var_table::Dict{Float64, String} # Map: height (m) => SMEAR API variable name (e.g. 16.8 => "HYY_META.T168")
end

# Pre-configured station maps for SMEAR-I and SMEAR-II tall towers
const SMEAR_I_VARRIO = SmearStationConfig(
    1, "SMEAR-I (Värriö)", 0.30, 15.0,
    Dict(2.0 => "VAR_META.T20", 8.0 => "VAR_META.T80", 15.0 => "VAR_META.T150")
)

const SMEAR_II_HYYTIALA = SmearStationConfig(
    2, "SMEAR-II (Hyytiälä)", 0.80, 125.0,
    Dict(
        4.2 => "HYY_META.T42", 8.4 => "HYY_META.T84", 16.8 => "HYY_META.T168",
        33.6 => "HYY_META.T336", 50.4 => "HYY_META.T504", 67.2 => "HYY_META.T672",
        125.0 => "HYY_META.T1250"
    )
)

"""
    fetch_smear_observations(config, t_start, t_end; timeout=30)

Queries the SMEAR REST API over a specified time window, handles HTTP/JSON parsing,
and pivots vertical profile streams into robust `StandardizedBLObservation` objects.
"""
function fetch_smear_observations(
    config::SmearStationConfig,
    t_start::DateTime,
    t_end::DateTime;
    timeout::Int = 30
)::Vector{StandardizedBLObservation}

    # 1. Construct endpoint URL
    table_variables = join(values(config.var_table), ",")
    url = "https://smear-backend.app.lab.helsinki.fi/v1/samedata?" *
          "tablevariable=$(table_variables)&" *
          "from=$(Dates.format(t_start, "yyyy-mm-ddTHH:MM:SS"))&" *
          "to=$(Dates.format(t_end, "yyyy-mm-ddTHH:MM:SS"))&" *
          "quality=ANY&aggregation=NONE"

    # 2. HTTP GET with timeout
    response = try
        HTTP.get(url; connect_timeout=timeout, readtimeout=timeout)
    catch e
        @error "SMEAR API request failed for $(config.name): $e"
        return StandardizedBLObservation[]
    end

    if response.status != 200
        @warn "SMEAR API returned HTTP status $(response.status)"
        return StandardizedBLObservation[]
    end

    # 3. Parse JSON payload
    json_data = JSON3.read(String(response.body))
    return parse_smear_json(json_data, config)
end

"""
    parse_smear_json(json_data, config)

Transforms raw JSON payload from SMEAR backend into standardized observations.
"""
function parse_smear_json(
    json_data,
    config::SmearStationConfig
)::Vector{StandardizedBLObservation}

    # Map variable back to height
    var_to_height = Dict(v => k for (k, v) in config.var_table)

    # Group observations by timestamp
    records_by_time = Dict{DateTime, Dict{Float64, Float64}}()

    for entry in json_data
        dt_str = string(entry[:samptime])
        dt = tryparse(DateTime, dt_str)
        dt === nothing && continue

        if !haskey(records_by_time, dt)
            records_by_time[dt] = Dict{Float64, Float64}()
        end

        for (var_name, h) in var_to_height
            if haskey(entry, Symbol(var_name))
                val = entry[Symbol(var_name)]
                if val !== nothing && isfinite(Float64(val))
                    records_by_time[dt][h] = Float64(val)
                end
            end
        end
    end

    observations = StandardizedBLObservation[]

    # Sort timestamps chronologically
    for dt in sort(collect(keys(records_by_time)))
        height_val_map = records_by_time[dt]
        isempty(height_val_map) && continue

        # Extract sorted heights and values
        sorted_heights = sort(collect(keys(height_val_map)))
        sorted_values = [height_val_map[h] for h in sorted_heights]

        n_valid = length(sorted_heights)
        # Degree-3 Chebyshev expansion requires at least N >= 4 levels
        is_robust = n_valid >= 4

        # SMEAR API baseline placeholders for surface fluxes if not queried concurrently
        ustar = NaN
        L_obukhov = NaN

        push!(observations, StandardizedBLObservation(
            dt,
            config.name,
            sorted_heights,
            sorted_values,
            ustar,
            L_obukhov,
            config.z0m,
            is_robust,
            n_valid
        ))
    end

    return observations
end

end # module SmearApiAdapter