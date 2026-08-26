#!/usr/bin/env julia
using LinearAlgebra
using Printf

Base.@kwdef struct SCMParams
    g::Float64 = 9.81
    theta0::Float64 = 290.0
    kappa::Float64 = 0.4
    epsilon::Float64 = 0.05
    l0::Float64 = 0.5
    B0_max::Float64 = 0.05
    delta_reg::Float64 = 0.01
    beta::Float64 = 5.0
    G0::Float64 = 0.5
    gamma_s::Float64 = 1.5
    r_s::Float64 = 0.1
    e_floor::Float64 = 1e-4
end

function solve_natural_cubic_spline(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    h = diff(x)
    d = ones(n)
    dl = zeros(n - 1)
    du = zeros(n - 1)
    b = zeros(n)

    for i in 2:(n-1)
        dl[i-1] = h[i-1] / 6.0
        d[i] = (h[i-1] + h[i]) / 3.0
        du[i] = h[i] / 6.0
        b[i] = (y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]
    end
    return Tridiagonal(dl, d, du) \ b
end

function eval_spline_deriv(x_grid::AbstractVector{<:Real}, y::AbstractVector{<:Real}, M::AbstractVector{<:Real}, x_val::Real)
    n = length(x_grid)
    idx = x_val <= x_grid[1] ? 1 : (x_val >= x_grid[end] ? n - 1 : max(1, min(n - 1, searchsortedlast(x_grid, x_val))))

    xi, xip1 = x_grid[idx], x_grid[idx+1]
    hi = xip1 - xi

    Ai = M[idx] / (6.0 * hi)
    Bi = M[idx+1] / (6.0 * hi)
    Ci = y[idx] / hi - M[idx] * hi / 6.0
    Di = y[idx+1] / hi - M[idx+1] * hi / 6.0

    return -3.0 * Ai * (xip1 - x_val)^2 + 3.0 * Bi * (x_val - xi)^2 - Ci + Di
end

# Interval-scanning bracket root finder
function find_coordinate_folds(z::AbstractVector{<:Real}, y_lnL::AbstractVector{<:Real}, M::AbstractVector{<:Real})
    g(x_val) = x_val * eval_spline_deriv(z, y_lnL, M, x_val) - 1.0
    folds = Float64[]
    g_vals = [g(zi) for zi in z]

    for i in 1:(length(z)-1)
        if g_vals[i] * g_vals[i+1] <= 0
            a, b = z[i], z[i+1]
            ga, gb = g_vals[i], g_vals[i+1]
            for _ in 1:60
                c = 0.5 * (a + b)
                gc = g(c)
                if abs(gc) < 1e-6 || (b - a) < 1e-6
                    push!(folds, c)
                    break
                end
                if ga * gc < 0
                    b, gb = c, gc
                else
                    a, ga = c, gc
                end
            end
        end
    end
    return folds
end

# Phase-space 2D Jacobian determinant (det J)
function compute_min_jacobian_det(e::Vector{Float64}, S::Vector{Float64}, N2::Vector{Float64}, p::SCMParams)
    min_det = Inf
    for i in 1:length(e)
        Ri_val = N2[i] / (S[i]^2)
        denom_B = (e[i]^2 + p.delta_reg^2)^2
        dB_de = p.B0_max * (2.0 * e[i] * p.delta_reg^2) / denom_B
        dDest_de = (3.0 * e[i]^2) / (p.l0 * (1.0 + p.beta * Ri_val))

        F_e = (p.l0 * S[i]^2 - dB_de - dDest_de) / p.epsilon
        dDest_dS = (-p.beta * e[i]^3 / (p.l0 * (1.0 + p.beta * Ri_val)^2)) * (-2.0 * N2[i] / (S[i]^3))
        F_S = (2.0 * p.l0 * e[i] * S[i] - dDest_dS) / p.epsilon

        G_e = -p.gamma_s * S[i]
        G_S = -p.gamma_s * e[i] - p.r_s

        det_J = F_e * G_S - F_S * G_e
        min_det = min(min_det, det_J)
    end
    return min_det
end

function main()
    p = SCMParams()
    N_z = 50
    z = [1.0 + 199.0 * (i / (N_z - 1))^1.5 for i in 0:(N_z-1)]

    t_max, dt = 43200.0, 0.1
    steps = Int(t_max / dt)
    n_sub = 10
    dt_sub = dt / n_sub

    e, S = fill(0.6, N_z), fill(1.2, N_z)
    N2, Ri, L_profile = zeros(N_z), zeros(N_z), zeros(N_z)
    log_steps = Int(3600.0 / dt)

    println("="^110)
    println("              GSPT EMERGENT DYNAMICAL FOLD & SADDLE-NODE SIMULATION")
    println("="^110)
    @printf("%-8s | %-12s | %-12s | %-12s | %-14s | %-14s\n",
        "Time (h)", "Surf Temp(K)", "Peak Shear", "Max TKE", "Min det(J)", "Emergent z_fold")
    println("-"^110)

    for step in 1:steps
        t = step * dt
        theta_surf = p.theta0 - 5.0 * (t / t_max)
        delta_theta = 2.0 + 8.0 * (t / t_max)
        h_inv = 30.0 + 50.0 * sqrt(t / t_max)

        fac = (p.g / p.theta0) * (delta_theta / h_inv)
        @. N2 = fac * exp(-z / h_inv)

        @. S = max(S + dt * (p.G0 - p.gamma_s * e * S - p.r_s * S), 1e-4)

        @. Ri = N2 / (S^2)
        for _ in 1:n_sub
            for i in 1:N_z
                de = (p.l0 * e[i] * S[i]^2 - p.B0_max * e[i]^2 / (e[i]^2 + p.delta_reg^2) -
                      e[i]^3 / (p.l0 * (1.0 + p.beta * Ri[i]))) / p.epsilon
                e[i] = max(e[i] + dt_sub * de, p.e_floor)
            end
        end

        if step % log_steps == 0
            for i in 1:N_z
                theta_v = theta_surf + delta_theta * (1.0 - exp(-z[i] / h_inv))
                dtheta_dz = (delta_theta / h_inv) * exp(-z[i] / h_inv)
                Km = p.l0 * sqrt(e[i])
                Kh = p.l0 * sqrt(e[i])

                u_star3 = (Km * S[i])^1.5
                flux_q = Kh * dtheta_dz
                L_profile[i] = (u_star3 * theta_v) / (p.kappa * p.g * max(flux_q, 1e-8))
            end

            y_lnL = log.(L_profile)
            M = solve_natural_cubic_spline(z, y_lnL)
            folds = find_coordinate_folds(z, y_lnL, M)
            min_detJ = compute_min_jacobian_det(e, S, N2, p)

            fold_str = isempty(folds) ? "None Detected" : @sprintf("%.4f m", folds[1])
            @printf("%8.1f | %12.4f | %12.4f | %12.4f | %14.4e | %-14s\n",
                t/3600.0, theta_surf, maximum(S), maximum(e), min_detJ, fold_str)
        end
    end
    println("="^110)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end