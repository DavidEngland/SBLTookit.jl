using LinearAlgebra

"""
    compute_eta(z_obs, y_obs, z_min, z_max, C, P; w_obs=nothing)

Compute continuous Hilbert-space modal coordinates

    η = C M a

from discrete tower observations.

`C` contains empirical mode coefficients in the Chebyshev basis.
`M` is the exact Chebyshev L² mass matrix for the first-kind weight.
`w_obs`, when supplied, contains diagonal observational weights.
"""
function compute_eta(
    z_obs::Vector{Float64},
    y_obs::Vector{Float64},
    z_min::Float64,
    z_max::Float64,
    C::Matrix{Float64},
    P::Int;
    w_obs::Union{Nothing,Vector{Float64}} = nothing,
)
    N = length(z_obs)

    @assert length(y_obs) == N
    @assert z_max > z_min
    @assert size(C, 2) == P + 1
    @assert size(C, 1) == 3

    # Physical → computational coordinate
    ξ = @. 2.0 * (z_obs - z_min) / (z_max - z_min) - 1.0

    # Guard against floating-point excursions outside [-1,1].
    ξ = clamp.(ξ, -1.0, 1.0)

    # Chebyshev Vandermonde/evaluation matrix
    B = Matrix{Float64}(undef, N, P + 1)

    θ = acos.(ξ)

    for k in 0:P
        @views B[:, k + 1] .= cos.(k .* θ)
    end

    # Weighted or unweighted least squares
    if w_obs === nothing
        a = B \ y_obs
    else
        @assert length(w_obs) == N
        @assert all(isfinite, w_obs)
        @assert all(>=(0.0), w_obs)

        sqrtW = sqrt.(w_obs)

        Bw = B .* sqrtW
        yw = y_obs .* sqrtW

        a = Bw \ yw
    end

    # Exact Chebyshev L² mass diagonal
    Mdiag = Vector{Float64}(undef, P + 1)
    Mdiag[1] = π

    if P >= 1
        Mdiag[2:end] .= π / 2
    end

    # Continuous Hilbert-space projection
    η = C * (Mdiag .* a)

    return η
end