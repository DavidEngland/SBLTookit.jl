module SmearAdapter

using DataFrames
using Dates
using ..CoreTypes

export profile_from_legacy, profile_to_legacy,
       standardized_from_legacy, profile_from_standardized,
       standardized_from_profile

function profile_from_legacy(
    datetime::DateTime,
    heights::Vector{Float64},
    values::Vector{Float64},
    zeta::Float64,
    ustar::Float64;
    reference_height::Float64=23.0,
    canopy_displacement::Float64=0.0,
)
    L = abs(zeta) > 1.0e-12 ? reference_height / zeta : NaN
    meta = ProfileMetadata(datetime, ustar, L, reference_height, canopy_displacement)
    return MeteorologicalProfile(meta, heights, values)
end

function profile_to_legacy(prof::MeteorologicalProfile, n_obs::Int)
    zeta = isfinite(prof.metadata.obukhov_length) && abs(prof.metadata.obukhov_length) > 1.0e-12 ?
           prof.metadata.reference_height / prof.metadata.obukhov_length : NaN
    return (
        datetime = prof.metadata.timestamp,
        heights = prof.heights,
        values = prof.values,
        n_obs = n_obs,
        zeta = zeta,
        ustar = prof.metadata.friction_velocity,
    )
end

function standardized_from_legacy(
    datetime::DateTime,
    heights::Vector{Float64},
    values::Vector{Float64},
    zeta::Float64,
    ustar::Float64;
    campaign::String="SMEAR",
    reference_height::Float64=23.0,
    z0m::Float64=NaN,
    robust_for_eta3::Union{Nothing,Bool}=nothing,
)
    L = abs(zeta) > 1.0e-12 ? reference_height / zeta : NaN
    n_valid = min(length(heights), length(values))
    eta3_ok = isnothing(robust_for_eta3) ? n_valid >= 3 : robust_for_eta3
    return StandardizedBLObservation(
        datetime,
        campaign,
        heights,
        values,
        ustar,
        L,
        z0m,
        eta3_ok,
        n_valid,
    )
end

function standardized_from_profile(
    prof::MeteorologicalProfile;
    campaign::String="SMEAR",
    z0m::Float64=NaN,
    robust_for_eta3::Union{Nothing,Bool}=nothing,
)
    n_valid = min(length(prof.heights), length(prof.values))
    eta3_ok = isnothing(robust_for_eta3) ? n_valid >= 3 : robust_for_eta3
    return StandardizedBLObservation(
        prof.metadata.timestamp,
        campaign,
        prof.heights,
        prof.values,
        prof.metadata.friction_velocity,
        prof.metadata.obukhov_length,
        z0m,
        eta3_ok,
        n_valid,
    )
end

function profile_from_standardized(
    obs::StandardizedBLObservation;
    reference_height::Float64=23.0,
    canopy_displacement::Float64=0.0,
)
    meta = ProfileMetadata(
        obs.datetime,
        obs.ustar,
        obs.L_obukhov,
        reference_height,
        canopy_displacement,
    )
    return MeteorologicalProfile(meta, obs.heights, obs.values)
end

end # module SmearAdapter
