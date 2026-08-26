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

function find_coordinate_fold(z::AbstractVector{<:Real}, y_lnL::AbstractVector{<:Real}, M::AbstractVector{<:Real})
    g(x_val) = x_val * eval_spline_deriv(z, y_lnL, M, x_val) - 1.0

    a_bound, b_bound = z[1], z[end]
    ga, gb = g(a_bound), g(b_bound)

    if ga * gb >= 0
        return NaN
    end

    for _ in 1:100
        c_val = 0.5 * (a_bound + b_bound)
        gc = g(c_val)

        if abs(gc) < 1e-5 || (b_bound - a_bound) < 1e-5
            return c_val
        end

        if ga * gc < 0
            b_bound, gb = c_val, gc
        else
            a_bound, ga = c_val, gc
        end
    end

    return 0.5 * (a_bound + b_bound)
end

function main()
    p = SCMParams()
    N_z = 50
    z = [1.0 + (199.0) * (i / N_z)^1.5 for i in 0:(N_z-1)]

    t_max = 43200.0
    dt = 0.1
    steps = Int(t_max / dt)
    n_sub = 5
    dt_sub = dt / n_sub

    # Pre-allocated state & work arrays
    e = fill(0.6, N_z)
    S = fill(1.2, N_z)
    N2 = zeros(N_z)
    Ri = zeros(N_z)
    L_profile = zeros(N_z)

    log_steps = Int(3600.0 / dt)

    println("="^115)
    println("                        GSPT 1D SINGLE-COLUMN MODEL COORDINATE FOLD SIMULATION")
    println("="^115)
    @printf("%-10s | %-12s | %-12s | %-12s | %-12s | %-12s | %-12s\n",
        "Time (h)", "Surf Temp(K)", "Inv Depth(m)", "Peak Shear", "Max TKE", "Est z_fold", "Exact z_fold")
    println("-"^115)

    for step in 1:steps
        t = step * dt

        theta_surf = p.theta0 - 5.0 * (t / t_max)
        delta_theta = 2.0 + 8.0 * (t / t_max)
        h_inv = 30.0 + 50.0 * sqrt(t / t_max)

        # 1. Update N2 in-place
        fac = (p.g / p.theta0) * (delta_theta / h_inv)
        @. N2 = fac * exp(-z / h_inv)

        # 2. Integrate slow shear in-place
        @. S = max(S + dt * (p.G0 - p.gamma_s * e * S - p.r_s * S), 1e-4)

        # 3. Fast TKE sub-stepping in-place (Zero allocations)
        @. Ri = N2 / (S^2)
        for _ in 1:n_sub
            for i in 1:N_z
                de = (p.l0 * e[i] * S[i]^2 - p.B0_max * e[i]^2 / (e[i]^2 + p.delta_reg^2) -
                      e[i]^3 / (p.l0 * (1.0 + p.beta * Ri[i]))) / p.epsilon
                e[i] = max(e[i] + dt_sub * de, p.e_floor)
            end
        end

        # 4. Hourly Logging
        if step % log_steps == 0
            a_param = 0.5 + 4.5 * (t / t_max)
            L0 = 15.0 + 10.0 * (t / t_max)
            @. L_profile = L0 * exp(a_param * z / 200.0)

            y_lnL = log.(L_profile)
            M = solve_natural_cubic_spline(z, y_lnL)
            z_fold_est = find_coordinate_fold(z, y_lnL, M)
            z_fold_exact = 200.0 / a_param

            z_fold_str = isnan(z_fold_est) ? "Out of Bound" : @sprintf("%.4f", z_fold_est)
            @printf("%10.1f | %12.4f | %12.4f | %12.4f | %12.4f | %12s | %12.4f\n",
                t/3600.0, theta_surf, h_inv, maximum(S), maximum(e), z_fold_str, z_fold_exact)
        end
    end
    println("="^115)
end

main()