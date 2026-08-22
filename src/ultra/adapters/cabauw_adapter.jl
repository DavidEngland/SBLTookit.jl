module CabauwAdapter

using DataFrames
using Dates
using ..CoreTypes

export CabauwTower, extract_temperature_profiles, extract_temperature_observations

struct CabauwTower <: AbstractObservationalTower
    z_nodes::Vector{Float64}
    reference_height::Float64
    displacement::Float64

    function CabauwTower()
        new([2.0, 10.0, 20.0, 40.0, 80.0, 140.0, 200.0], 20.0, 0.0)
    end
end

function extract_temperature_profiles(df::DataFrame, tower::CabauwTower=CabauwTower())
    profiles = MeteorologicalProfile[]
    cabauw_t_cols = ["TA_2m", "TA_10m", "TA_20m", "TA_40m", "TA_80m", "TA_140m", "TA_200m"]

    missing_cols = filter(c -> !(c in names(df)), cabauw_t_cols)
    isempty(missing_cols) || error("Cabauw input missing columns: $(missing_cols)")

    for row in eachrow(df)
        raw_vals = Any[row[c] for c in cabauw_t_cols]
        valid_mask = falses(length(raw_vals))
        vals = Float64[]
        z = Float64[]

        for i in eachindex(raw_vals)
            v = raw_vals[i]
            if !ismissing(v)
                fv = tryparse(Float64, string(v))
                if !isnothing(fv) && isfinite(fv)
                    valid_mask[i] = true
                    push!(vals, fv)
                    push!(z, tower.z_nodes[i])
                end
            end
        end

        length(vals) >= 3 || continue

        ustar = row_hasproperty(row, :ustar) ? to_float(row.ustar, NaN) : NaN
        L = row_hasproperty(row, :L_obukhov) ? to_float(row.L_obukhov, NaN) : NaN
        timestamp = row_hasproperty(row, :datetime) ? DateTime(row.datetime) : continue

        meta = ProfileMetadata(timestamp, ustar, L, tower.reference_height, tower.displacement)
        push!(profiles, MeteorologicalProfile(meta, z, vals))
    end

    return profiles
end

function extract_temperature_observations(
    df::DataFrame,
    tower::CabauwTower=CabauwTower();
    campaign::String="GABLS3",
    z0m::Float64=0.15,
)
    profiles = extract_temperature_profiles(df, tower)
    observations = StandardizedBLObservation[]

    for prof in profiles
        n_valid = min(length(prof.heights), length(prof.values))
        push!(observations, StandardizedBLObservation(
            prof.metadata.timestamp,
            campaign,
            prof.heights,
            prof.values,
            prof.metadata.friction_velocity,
            prof.metadata.obukhov_length,
            z0m,
            n_valid >= 3,
            n_valid,
        ))
    end

    return observations
end

@inline function row_hasproperty(row, name::Symbol)
    return name in propertynames(row)
end

@inline function to_float(v, default::Float64)
    if v isa Number
        return Float64(v)
    elseif ismissing(v)
        return default
    end
    parsed = tryparse(Float64, string(v))
    return isnothing(parsed) ? default : parsed
end

end # module CabauwAdapter
