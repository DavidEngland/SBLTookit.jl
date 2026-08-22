module ICOSAdapter

using Dates
using ..CoreTypes: MeteorologicalProfile, ProfileMetadata, StandardizedBLObservation

export upscale_sparse_icos_profile, upscale_sparse_icos_observation

@inline function _psi_m(zeta::Float64)
    if zeta <= 0.0
        x = (1.0 - 16.0 * zeta)^0.25
        return 2.0 * log((1.0 + x) / 2.0) + log((1.0 + x^2) / 2.0) - 2.0 * atan(x) + pi / 2.0
    end
    return -5.0 * zeta
end

@inline function _psi_h(zeta::Float64)
    if zeta <= 0.0
        y = sqrt(1.0 - 16.0 * zeta)
        return 2.0 * log((1.0 + y) / 2.0)
    end
    return -5.0 * zeta
end

@inline function _stability_coordinate(z::Float64, L::Float64, stability_correction::Symbol)
    chi = log(z)
    if stability_correction == :none
        return chi
    elseif !isfinite(L) || abs(L) <= 1.0e-9
        # Undefined or near-singular Obukhov scale: fall back to neutral log coordinate.
        return chi
    end

    zeta = z / L
    if stability_correction == :psi_m
        return chi - _psi_m(zeta)
    elseif stability_correction == :psi_h
        return chi - _psi_h(zeta)
    end

    error("Unsupported stability_correction=$(stability_correction). Use :none, :psi_m, or :psi_h")
end

"""
    upscale_sparse_icos_profile(dt, z_low, z_high, v_low, v_high, ustar, L;
                                target_levels=[2.0, 5.0, 10.0, 20.0, 40.0],
                                stability_correction=:none)

Interpolate sparse two-level tower observations into a denser virtual profile
using log-linear surface-layer scaling. Optional Monin-Obukhov stability
corrections can be enabled with `stability_correction=:psi_m` or `:psi_h`.
"""
function upscale_sparse_icos_profile(
    dt::DateTime,
    z_low::Float64,
    z_high::Float64,
    v_low::Float64,
    v_high::Float64,
    ustar::Float64,
    L::Float64;
    target_levels::Vector{Float64}=[2.0, 5.0, 10.0, 20.0, 40.0],
    stability_correction::Symbol=:none,
)
    z_low > 0.0 || error("z_low must be above ground")
    z_high > z_low || error("z_high must be greater than z_low")

    χ_low = _stability_coordinate(z_low, L, stability_correction)
    χ_high = _stability_coordinate(z_high, L, stability_correction)
    denom = χ_high - χ_low
    abs(denom) > 1e-6 || error("Levels are too close to resolve a log slope")

    slope = (v_high - v_low) / denom

    virtual_vals = Float64[]
    for z_t in target_levels
        z_t > 0.0 || error("target_levels must be above ground")
        χ_t = _stability_coordinate(z_t, L, stability_correction)
        v_t = v_low + slope * (χ_t - χ_low)
        push!(virtual_vals, v_t)
    end

    meta = ProfileMetadata(dt, ustar, L, z_high, 0.0)
    return MeteorologicalProfile(meta, target_levels, virtual_vals)
end

"""
    upscale_sparse_icos_observation(...)

Generate an upscaled ICOS profile and return a standardized campaign-normalized
observation payload for unified ingestion.
"""
function upscale_sparse_icos_observation(
    dt::DateTime,
    z_low::Float64,
    z_high::Float64,
    v_low::Float64,
    v_high::Float64,
    ustar::Float64,
    L::Float64;
    campaign::String="ICOS-Upscaled",
    z0m::Float64=NaN,
    target_levels::Vector{Float64}=[2.0, 5.0, 10.0, 20.0, 40.0],
    stability_correction::Symbol=:none,
)
    prof = upscale_sparse_icos_profile(
        dt,
        z_low,
        z_high,
        v_low,
        v_high,
        ustar,
        L;
        target_levels,
        stability_correction,
    )

    n_valid = min(length(prof.heights), length(prof.values))
    return StandardizedBLObservation(
        prof.metadata.timestamp,
        campaign,
        prof.heights,
        prof.values,
        prof.metadata.friction_velocity,
        prof.metadata.obukhov_length,
        z0m,
        n_valid >= 3,
        n_valid,
    )
end

end # module ICOSAdapter
