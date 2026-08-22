module CoreTypes

using Dates

export AbstractObservationalTower, ProfileMetadata, MeteorologicalProfile, StandardizedBLObservation

abstract type AbstractObservationalTower end

struct ProfileMetadata
    timestamp::DateTime
    friction_velocity::Float64
    obukhov_length::Float64
    reference_height::Float64
    canopy_displacement::Float64
end

struct MeteorologicalProfile
    metadata::ProfileMetadata
    heights::Vector{Float64}
    values::Vector{Float64}

    function MeteorologicalProfile(meta::ProfileMetadata, heights::Vector{Float64}, values::Vector{Float64})
        @assert length(heights) == length(values) "Height and value vectors must have equal length."
        @assert issorted(heights) "Tower heights must be strictly monotonically increasing."
        return new(meta, heights, values)
    end
end

struct StandardizedBLObservation
    datetime::DateTime
    campaign::String
    heights::Vector{Float64}
    values::Vector{Float64}
    ustar::Float64
    L_obukhov::Float64
    z0m::Float64
    robust_for_eta3::Bool
    n_valid_levels::Int

    function StandardizedBLObservation(
        datetime::DateTime, campaign::String, heights::Vector{Float64}, values::Vector{Float64},
        ustar::Float64, L_obukhov::Float64, z0m::Float64, min_levels::Int=4
    )
        @assert length(heights) == length(values) "Heights and values length mismatch."
        @assert issorted(heights) "Heights must be strictly increasing."

        n_valid = length(heights)
        # Verify valid physics and sufficient sensor count for P >= 3 degree fits
        is_robust = (n_valid >= min_levels) && isfinite(ustar) && (ustar > 0.0) && isfinite(L_obukhov)

        return new(datetime, campaign, heights, values, ustar, L_obukhov, z0m, is_robust, n_valid)
    end
end

end # module CoreTypes