# In src/ultra/adapters/smear.jl:
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

    # Filter out invalid entries before counting
    valid_mask = .!(isnan.(heights) .| isnan.(values))
    n_valid = count(valid_mask)

    # Require N >= 4 for robust degree-3 polynomial fits
    eta3_ok = isnothing(robust_for_eta3) ? (n_valid >= 4 && isfinite(ustar) && ustar > 0.0) : robust_for_eta3

    return StandardizedBLObservation(
        datetime,
        campaign,
        heights[valid_mask],
        values[valid_mask],
        ustar,
        L,
        z0m,
        eta3_ok,
        n_valid,
    )
end