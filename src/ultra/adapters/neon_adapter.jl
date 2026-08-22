module NEONAdapter

using Dates
using DataFrames
using ..CoreTypes: MeteorologicalProfile, ProfileMetadata, StandardizedBLObservation

export extract_neon_profile, extract_neon_observations

"""
    extract_neon_profile(df, prefix, heights; reference_height=10.0)

Extract timestamped profiles from unstacked NEON-like columns named
`<prefix>_z1`, `<prefix>_z2`, ..., `<prefix>_zN` where N is `length(heights)`.
Rows missing any target level are skipped.
"""
function extract_neon_profile(
    df::DataFrame,
    prefix::String,
    heights::Vector{Float64};
    reference_height::Float64=10.0,
)
    n_levels = length(heights)
    n_levels >= 2 || error("Need at least two heights for a profile")

    col_names = ["$(prefix)_z$(i)" for i in 1:n_levels]
    available = filter(c -> c in names(df), col_names)
    length(available) == n_levels || error("Missing expected NEON profile columns. Found: $(available)")

    datetime_col = :datetime in names(df) ? :datetime : (:timestamp in names(df) ? :timestamp : nothing)
    datetime_col === nothing && error("Expected a datetime/timestamp column in NEON input")

    profiles = MeteorologicalProfile[]
    for row in eachrow(df)
        dt_raw = row[datetime_col]
        dt = dt_raw isa DateTime ? dt_raw : tryparse(DateTime, string(dt_raw))
        dt === nothing && continue

        vals = Float64[]
        valid = true
        for cname in col_names
            v = row[Symbol(cname)]
            if ismissing(v)
                valid = false
                break
            end
            fv = v isa Number ? Float64(v) : something(tryparse(Float64, string(v)), NaN)
            if !isfinite(fv)
                valid = false
                break
            end
            push!(vals, fv)
        end
        valid || continue

        ustar = _row_float(row, [:ustar, :u_star, :friction_velocity], NaN)
        L = _row_float(row, [:L, :L_obukhov, :mo_length], NaN)

        meta = ProfileMetadata(dt, ustar, L, reference_height, 0.0)
        push!(profiles, MeteorologicalProfile(meta, heights, vals))
    end

    return profiles
end

"""
    extract_neon_observations(df, prefix, heights; campaign="NEON", z0m=NaN,
                              min_levels=2, reference_height=10.0)

Extract campaign-normalized observations from NEON-like unstacked columns
`<prefix>_z1`, `<prefix>_z2`, ..., `<prefix>_zN`. Unlike `extract_neon_profile`,
this routine allows partially populated levels and records per-row level coverage.
"""
function extract_neon_observations(
    df::DataFrame,
    prefix::String,
    heights::Vector{Float64};
    campaign::String="NEON",
    z0m::Float64=NaN,
    min_levels::Int=2,
    reference_height::Float64=10.0,
)
    n_levels = length(heights)
    n_levels >= 2 || error("Need at least two heights for a profile")

    col_names = ["$(prefix)_z$(i)" for i in 1:n_levels]
    available = filter(c -> c in names(df), col_names)
    length(available) == n_levels || error("Missing expected NEON profile columns. Found: $(available)")

    datetime_col = :datetime in names(df) ? :datetime : (:timestamp in names(df) ? :timestamp : nothing)
    datetime_col === nothing && error("Expected a datetime/timestamp column in NEON input")

    observations = StandardizedBLObservation[]
    for row in eachrow(df)
        dt_raw = row[datetime_col]
        dt = dt_raw isa DateTime ? dt_raw : tryparse(DateTime, string(dt_raw))
        dt === nothing && continue

        vals = Float64[]
        z = Float64[]
        for (i, cname) in enumerate(col_names)
            v = row[Symbol(cname)]
            ismissing(v) && continue
            fv = v isa Number ? Float64(v) : something(tryparse(Float64, string(v)), NaN)
            isfinite(fv) || continue
            push!(z, heights[i])
            push!(vals, fv)
        end

        n_valid = min(length(z), length(vals))
        n_valid >= min_levels || continue

        ustar = _row_float(row, [:ustar, :u_star, :friction_velocity], NaN)
        L = _row_float(row, [:L, :L_obukhov, :mo_length], NaN)

        # Keep metadata parity with profile path; reference height is reserved for
        # downstream reconstruction when needed.
        _ = reference_height
        push!(observations, StandardizedBLObservation(
            dt,
            campaign,
            z,
            vals,
            ustar,
            L,
            z0m,
            n_valid >= 3,
            n_valid,
        ))
    end

    return observations
end

@inline function _row_float(row, candidates::Vector{Symbol}, default::Float64)
    row_names = propertynames(row)
    for col in candidates
        if col in row_names
            value = row[col]
            if value isa Number
                return Float64(value)
            elseif ismissing(value)
                continue
            end
            parsed = tryparse(Float64, string(value))
            parsed === nothing || return parsed
        end
    end
    return default
end

end # module NEONAdapter
