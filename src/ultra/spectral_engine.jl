module SpectralEngine

using LinearAlgebra
using ..CoreTypes

export StaticProjector, chebyshev_fingerprint

struct StaticProjector
    heights::Vector{Float64}
    n_coeffs::Int
    height_mapping::Symbol
    pseudo_inverse::Matrix{Float64}
end

const _PROJECTOR_CACHE = Dict{Tuple{Tuple{Vararg{Float64}}, Int, Symbol}, StaticProjector}()

function _map_to_chebyshev(heights::Vector{Float64}, height_mapping::Symbol)
    h_min, h_max = extrema(heights)
    h_max > h_min || error("Heights must span a non-zero interval")

    if height_mapping == :linear
        return @. 2.0 * (heights - h_min) / (h_max - h_min) - 1.0
    elseif height_mapping == :log
        all(>(0.0), heights) || error("Log-mapping requires strictly positive heights")
        return @. 2.0 * (log(heights) - log(h_min)) / (log(h_max) - log(h_min)) - 1.0
    end
    error("Unknown height mapping: $(height_mapping)")
end

function _get_projector(heights::Vector{Float64}, n_coeffs::Int, height_mapping::Symbol)
    key = (Tuple(heights), n_coeffs, height_mapping)
    return get!(_PROJECTOR_CACHE, key) do
        x = _map_to_chebyshev(heights, height_mapping)
        deg = min(n_coeffs - 1, length(heights) - 1)
        V = zeros(Float64, length(x), deg + 1)
        for (j, xj) in enumerate(x)
            theta = acos(clamp(xj, -1.0, 1.0))
            for k in 0:deg
                V[j, k + 1] = cos(k * theta)
            end
        end
        pinv_v = pinv(V)
        padded = zeros(Float64, n_coeffs, length(heights))
        padded[1:size(pinv_v, 1), :] .= pinv_v
        StaticProjector(copy(heights), n_coeffs, height_mapping, padded)
    end
end

function chebyshev_fingerprint(
    prof::MeteorologicalProfile;
    n_coeffs::Int=4,
    height_mapping::Symbol=:log,
)
    adjusted_heights = prof.heights .- prof.metadata.canopy_displacement
    length(adjusted_heights) == length(prof.values) || error("heights and values must be same length")
    length(adjusted_heights) >= 2 || error("Need >= 2 points for Chebyshev transform")

    projector = _get_projector(adjusted_heights, n_coeffs, height_mapping)
    result = zeros(Float64, n_coeffs)
    mul!(result, projector.pseudo_inverse, prof.values)
    return result
end

end # module SpectralEngine
