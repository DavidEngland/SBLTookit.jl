#!/usr/bin/env julia
# Computes mesh-invariant pattern coordinates eta from observational tower data using Chebyshev polynomial fitting and Hilbert space projection.
using LinearAlgebra

"""
    compute_eta(z_obs, y_obs, z_min, z_max, C, P)

Computes mesh-invariant pattern coordinates eta from observation points.
"""
function compute_eta(
    z_obs::Vector{Float64},
    y_obs::Vector{Float64},
    z_min::Float64,
    z_max::Float64,
    C::Matrix{Float64},
    P::Int
)
    # Step 1: Normalize heights to [-1, 1]
    xi = @. 2.0 * (z_obs - z_min) / (z_max - z_min) - 1.0

    # Step 2: Build evaluation matrix B
    N = length(z_obs)
    B = zeros(N, P + 1)
    for k in 0:P
        B[:, k + 1] = cos.(k .* acos.(xi))
    end

    # Fit polynomial coefficients via least-squares
    a = B \ y_obs

    # Step 3: Construct diagonal Mass Matrix M
    M_diag = [k == 0 ? π : π / 2.0 for k in 0:P]
    M = Diagonal(M_diag)

    # Steps 4 & 5: Project onto empirical mode matrix C
    eta = C * M * a
    return eta
end