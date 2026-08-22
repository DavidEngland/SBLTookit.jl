module ShebaAdapter

using DataFrames
using Dates
using ..CoreTypes

export ShebaTower, extract_sheba_profiles, extract_sheba_observations

struct ShebaTower <: AbstractObservationalTower
    z_lo::Float64
    z_hi::Float64
    reference_height::Float64
    displacement::Float64

    function ShebaTower(; z_lo::Float64=2.5, z_hi::Float64=10.0)
        new(z_lo, z_hi, z_hi, 0.0)
    end
end

function extract_sheba_profiles(
    df::DataFrame;
    value_pair::Symbol=:temperature,
    z_lo::Float64=2.5,
    z_hi::Float64=10.0,
)::Vector{MeteorologicalProfile}
    col_pairs = Dict(
        :temperature => (:T_lo, :T_hi),
        :wind => (:ws_lo, :ws_hi),
        :humidity => (:q_lo, :q_hi),
    )

    haskey(col_pairs, value_pair) || error("value_pair must be one of: $(collect(keys(col_pairs)))")

    lo_col, hi_col = col_pairs[value_pair]
    missing_cols = filter(c -> !(c in names(df)), [lo_col, hi_col])
    isempty(missing_cols) || error("Input missing columns: $(missing_cols)")

    tower = ShebaTower(; z_lo, z_hi)
    profiles = MeteorologicalProfile[]

    for row in eachrow(df)
        lo_val = _to_float(row[lo_col], NaN)
        hi_val = _to_float(row[hi_col], NaN)
        isfinite(lo_val) && isfinite(hi_val) || continue

        timestamp = _row_datetime(row)
        timestamp === nothing && continue

        ustar = _row_float(row, :ustar, NaN)
        L = _row_float(row, :L_obukhov, NaN)
        meta = ProfileMetadata(timestamp, ustar, L, tower.reference_height, tower.displacement)
        push!(profiles, MeteorologicalProfile(meta, [tower.z_lo, tower.z_hi], [lo_val, hi_val]))
    end

    return profiles
end

function extract_sheba_observations(
    df::DataFrame;
    campaign::String="SHEBA",
    z0m::Float64=0.0001,
    value_pair::Symbol=:temperature,
    z_lo::Float64=2.5,
    z_hi::Float64=10.0,
)::Vector{StandardizedBLObservation}
    profiles = extract_sheba_profiles(df; value_pair, z_lo, z_hi)
    observations = StandardizedBLObservation[]

    for prof in profiles
        n_valid = min(count(isfinite, prof.heights), count(isfinite, prof.values))
        push!(observations, StandardizedBLObservation(
            prof.metadata.timestamp,
            campaign,
            prof.heights,
            prof.values,
            prof.metadata.friction_velocity,
            prof.metadata.obukhov_length,
            z0m,
            false,
            n_valid,
        ))
    end

    return observations
end

@inline function _row_datetime(row)
    for col in (:datetime, :timestamp)
        if col in propertynames(row)
            val = row[col]
            val isa DateTime && return val
            parsed = tryparse(DateTime, string(val))
            return parsed
        end
    end
    return nothing
end

@inline function _row_float(row, col::Symbol, default::Float64)
    col in propertynames(row) || return default
    return _to_float(row[col], default)
end

@inline function _to_float(v, default::Float64)
    if v isa Number
        f = Float64(v)
        return isfinite(f) ? f : default
    elseif ismissing(v)
        return default
    end
    parsed = tryparse(Float64, string(v))
    return isnothing(parsed) ? default : (isfinite(parsed) ? parsed : default)
end

end # module ShebaAdapter
