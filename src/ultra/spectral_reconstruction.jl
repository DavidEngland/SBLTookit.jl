module SpectralReconstruction

using LinearAlgebra

export reconstruct_from_chebyshev

"""
    reconstruct_from_chebyshev(
        coeffs::Vector{Float64},
        target_heights::Vector{Float64},
        height_bounds::Tuple{Float64, Float64};
        height_mapping::Symbol=:log,
        canopy_displacement::Float64=0.0
    ) -> Vector{Float64}

Synthesizes vertical atmospheric field values on `target_heights` from Chebyshev modes.
`height_bounds` must match the `(z_min, z_max)` interval used in forward projection.
"""
function reconstruct_from_chebyshev(
    coeffs::Vector{Float64},
    target_heights::Vector{Float64},
    height_bounds::Tuple{Float64, Float64};
    height_mapping::Symbol=:log,
    canopy_displacement::Float64=0.0
)
    # Apply zero-plane canopy displacement shift
    adj_targets = target_heights .- canopy_displacement
    h_min, h_max = height_bounds .- canopy_displacement

    h_max > h_min || error("Height bounds must span a positive non-zero interval")
    all(>(0.0), adj_targets) || error("Target heights after displacement must be strictly positive")

    # Map physical target heights to Chebyshev domain [-1.0, 1.0]
    if height_mapping == :linear
        x = @. 2.0 * (adj_targets - h_min) / (h_max - h_min) - 1.0
    elseif height_mapping == :log
        x = @. 2.0 * (log(adj_targets) - log(h_min)) / (log(h_max) - log(h_min)) - 1.0
    else
        error("Unknown height mapping: $(height_mapping)")
    end

    # Clamp domain to guard against numerical overflow past target boundaries
    clamp!(x, -1.0, 1.0)

    # Build evaluation Vandermonde matrix on grid targets
    n_targets = length(target_heights)
    n_coeffs = length(coeffs)
    V_target = zeros(Float64, n_targets, n_coeffs)

    for (i, xi) in enumerate(x)
        theta = acos(xi)
        for k in 0:(n_coeffs - 1)
            V_target[i, k + 1] = cos(k * theta)
        end
    end

    # Synthesize vertical field profile
    return V_target * coeffs
end

end # module SpectralReconstruction