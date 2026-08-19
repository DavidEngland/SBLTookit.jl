#!/usr/bin/env julia
# SBLTookit.jl - Transformations for SBL data, including Planar Fit for 3D sonic anemometer winds.
# src/Transformations.jl
using LinearAlgebra

"""
Planar Fit transformation for raw 3D sonic anemometer winds over tilted/snowy terrain.
Solves w_bar = b0 + b1*u_bar + b2*v_bar to construct the rotation matrix P.
"""
function apply_planar_fit(u::Vector{Float64}, v::Vector{Float64}, w::Vector{Float64})
    N = length(u)
    # Regression matrix X = [1, u, v]
    X = hcat(ones(N), u, v)
    b = X \ w  # Solve b = [b0, b1, b2]

    b0, b1, b2 = b[1], b[2], b[3]

    # Calculate pitch (alpha) and roll (beta) angles
    p = tan(atan(b1))
    q = tan(atan(b2))
    cos_gamma = 1.0 / sqrt(1.0 + b1^2 + b2^2)

    sin_alpha = -b1 * cos_gamma
    cos_alpha = sqrt(1.0 - sin_alpha^2)

    sin_beta = b2 / sqrt(1.0 + b2^2)
    cos_beta = sqrt(1.0 - sin_beta^2)

    # Rotation matrix
    R = [
         cos_alpha              0.0         -sin_alpha;
        -sin_beta*sin_alpha   cos_beta     -sin_beta*cos_alpha;
         cos_beta*sin_alpha   sin_beta      cos_beta*cos_alpha
    ]

    # Rotate velocity vectors
    V_raw = hcat(u, v, w)'
    V_rot = R * V_raw

    u_pf = V_rot[1, :]
    v_pf = V_rot[2, :]
    w_pf = V_rot[3, :] .- b0 # Subtract residual bias

    return u_pf, v_pf, w_pf
end # module Transformations