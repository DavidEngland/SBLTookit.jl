module BLLASTAdapter

using Dates
using DataFrames
using ..CoreTypes

export standardized_from_bllast, profile_from_bllast,
       bllast_df_to_standardized, BLLASTSiteMetadata

"""
    BLLASTSiteMetadata

Default physical parameters for BLLAST campaign sites (e.g., 60m tower, 120m site).
"""
struct BLLASTSiteMetadata
    site_id::String
    latitude::Float64
    longitude::Float64
    reference_height::Float64
    z0m::Float64
    canopy_displacement::Float64
end

# Pre-configured metadata for the main BLLAST 60m instrumented tower (Lannemezan)
const BLLAST_60M_SITE = BLLASTSiteMetadata("BLLAST_60M", 43.1238, 0.3621, 10.0, 0.15, 0.0)

"""
    standardized_from_bllast(datetime, raw_heights, raw_values, ustar, L_obukhov; kwargs...)

Cleans raw BLLAST tower profiles by filtering NaNs, sorting heights, and
mapping to `StandardizedBLObservation`.
"""
function standardized_from_bllast(
    datetime::DateTime,
    raw_heights::Vector{Float64},
    raw_values::Vector{Float64},
    ustar::Float64,
    L_obukhov::Float64;
    campaign::String = "BLLAST",
    site_meta::BLLASTSiteMetadata = BLLAST_60M_SITE,
    min_levels::Int = 4,
)
    # 1. Filter out NaNs, missing readings, or non-physical values
    valid_mask = .!(isnan.(raw_heights) .| isnan.(raw_values))
    heights = raw_heights[valid_mask]
    values = raw_values[valid_mask]

    # 2. Enforce monotonic height ordering required for p-FEM expansion
    p = sortperm(heights)
    sorted_heights = heights[p]
    sorted_values = values[p]

    n_valid = length(sorted_heights)

    # 3. Assess robustness (degree P=3 Chebyshev fits require N >= 4 levels)
    is_robust = (n_valid >= min_levels) && isfinite(ustar) && (ustar > 0.0) && isfinite(L_obukhov)

    return StandardizedBLObservation(
        datetime,
        campaign,
        sorted_heights,
        sorted_values,
        ustar,
        L_obukhov,
        site_meta.z0m,
        is_robust,
        n_valid
    )
end

"""
    bllast_df_to_standardized(df_row, height_cols; kwargs...)

Ingests a single row from a BLLAST processed DataFrame/CSV payload into `StandardizedBLObservation`.
`height_cols` maps target height values to DataFrame column symbols (e.g., `Dict(10.0 => :temp_10m, 60.0 => :temp_60m)`).
"""
function bllast_df_to_standardized(
    df_row::DataFrameRow,
    height_cols::Dict{Float64, Symbol};
    datetime_col::Symbol = :datetime,
    ustar_col::Symbol = :ustar,
    L_col::Symbol = :L_obukhov,
    site_meta::BLLASTSiteMetadata = BLLAST_60M_SITE
)
    raw_heights = Float64[]
    raw_values = Float64[]

    for (z, col_name) in height_cols
        val = Float64(df_row[col_name])
        push!(raw_heights, z)
        push!(raw_values, val)
    end

    dt = df_row[datetime_col]
    ustar = Float64(df_row[ustar_col])
    L = Float64(df_row[L_col])

    return standardized_from_bllast(dt, raw_heights, raw_values, ustar, L; site_meta=site_meta)
end

"""
    profile_from_bllast(obs; site_meta)

Converts a standardized BLLAST payload back to a `MeteorologicalProfile`.
"""
function profile_from_bllast(
    obs::StandardizedBLObservation;
    site_meta::BLLASTSiteMetadata = BLLAST_60M_SITE
)
    meta = ProfileMetadata(
        obs.datetime,
        obs.ustar,
        obs.L_obukhov,
        site_meta.reference_height,
        site_meta.canopy_displacement
    )
    return MeteorologicalProfile(meta, obs.heights, obs.values)
end

end # module BLLASTAdapter