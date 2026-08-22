module UnifiedBLIngestion

using DataFrames
using Dates
using ..CoreTypes: StandardizedBLObservation, MeteorologicalProfile, ProfileMetadata
using ..CabauwAdapter
using ..AmeriFluxAdapter
using ..ICOSAdapter
using ..NEONAdapter
using ..SpectralEngine

export detect_network_format, ingest_boundary_layer_data, batch_spectral_fingerprints

"""
    detect_network_format(df::DataFrame) -> Symbol

Inspects DataFrame schema to auto-classify the boundary layer network data source.
"""
function detect_network_format(df::DataFrame)::Symbol
    cols = string.(names(df))

    if "TA_2m" in cols && "TA_200m" in cols
        return :cabauw
    elseif "TIMESTAMP_START" in cols || any(c -> startswith(c, "TA_1_") || startswith(c, "TA_2_"), cols)
        return :ameriflux
    elseif any(c -> occursin(r"_z\d+$", c), cols)
        return :neon
    elseif "T_10m" in cols && "T_50m" in cols
        return :icos
    else
        return :unknown
    end
end

"""
    ingest_boundary_layer_data(df::DataFrame; kwargs...) -> Vector{StandardizedBLObservation}

Universal entrypoint: ingests raw observational DataFrames and normalizes them into
`StandardizedBLObservation` objects.
"""
function ingest_boundary_layer_data(
    df::DataFrame;
    site_id::String="",
    prefix::String="temp",
    heights::Vector{Float64}=Float64[],
    icos_target_levels::Vector{Float64}=[2.0, 5.0, 10.0, 20.0, 40.0],
    campaign::String="UNIFIED",
)::Vector{StandardizedBLObservation}
    format = detect_network_format(df)

    if format == :cabauw
        return extract_temperature_observations(df)

    elseif format == :neon
        isempty(heights) && error("NEON network format requires explicit `heights` argument.")
        return extract_neon_observations(df, prefix, heights; campaign=campaign)

    elseif format == :icos
        observations = StandardizedBLObservation[]
        for row in eachrow(df)
            dt = row.datetime isa DateTime ? row.datetime : DateTime(string(row.datetime))
            obs = upscale_sparse_icos_observation(
                dt, 10.0, 50.0, Float64(row.T_10m), Float64(row.T_50m),
                Float64(get(row, :ustar, 0.2)), Float64(get(row, :L, -50.0));
                campaign=campaign, target_levels=icos_target_levels, stability_correction=:psi_h
            )
            push!(observations, obs)
        end
        return observations

    elseif format == :ameriflux
        site_id != "" || error("AmeriFlux format requires a valid `site_id` parameter.")
        # Assumes file-path based processing for AmeriFlux CSVs
        error("Use `extract_ameriflux_observations(site_id, csv_path)` for AmeriFlux disk ingestion.")

    else
        error("Unable to auto-detect network schema from DataFrame columns.")
    end
end

"""
    batch_spectral_fingerprints(observations::Vector{StandardizedBLObservation}; n_coeffs=4)

Converts a list of standardized observations into Chebyshev spectral fingerprints.
"""
function batch_spectral_fingerprints(
    observations::Vector{StandardizedBLObservation};
    n_coeffs::Int=4
)
    matrix_out = zeros(Float64, length(observations), n_coeffs)

    for (i, obs) in enumerate(observations)
        meta = ProfileMetadata(obs.datetime, obs.ustar, obs.L_obukhov, obs.heights[end], 0.0)
        prof = MeteorologicalProfile(meta, obs.heights, obs.values)
        matrix_out[i, :] .= chebyshev_fingerprint(prof; n_coeffs=n_coeffs, height_mapping=:log)
    end

    return matrix_out
end

end # module UnifiedBLIngestion