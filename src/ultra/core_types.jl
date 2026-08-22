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
end

"""
    StandardizedBLObservation

Campaign-agnostic intermediate ingestion payload used to normalize profile adapters
before projection/spectral diagnostics.
"""
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
end

end # module CoreTypes
