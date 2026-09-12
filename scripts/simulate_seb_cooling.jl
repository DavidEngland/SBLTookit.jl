#!/usr/bin/env julia
# =============================================================================
# 1D SEB TIME-STEPPING SIMULATION: SPURIOUS QUENCHING vs C-Z0HR REGULARIZATION
# =============================================================================

using Printf
using LinearAlgebra

struct SEBModelConfig{T<:AbstractFloat}
    # Atmospheric Grid & Time Parameters
    Nz::Int           # Number of vertical grid levels (default: 50)
    z_max::T          # Top of column [m] (default: 200.0)
    dt::T             # Timestep [s] (default: 60.0)
    t_end::T          # Total integration time [s] (default: 12 * 3600.0)

    # Surface Energy Budget Parameters
    Cs::T             # Volumetric heat capacity of soil slab [J/m²/K] (2.0e5)
    Ks::T             # Deep soil conductivity coefficient [W/m²/K] (2.0)
    T_deep::T         # Deep soil temperature [K] (285.0)
    Rn0::T            # Net surface radiative cooling [W/m²] (-50.0)
    rho_cp::T         # Air heat capacity constant [J/m³/K] (1200.0)

    # Closure Parameters
    Ri_c::T           # Critical Richardson number (0.20)
    alpha::T          # C-Z0HR damping control parameter (2.0)
    epsilon_c::T      # Curvature softening cap (0.05)
    eps_hyper::T      # C^∞ hyperbolic smoothing width (1e-3)
    l_m::T            # Mixing length scale [m] (15.0)
    S_eff::T          # Constant ambient wind shear [s⁻¹] (0.015)
end

function SEBModelConfig(;
    Nz::Int = 50,
    z_max::Float64 = 200.0,
    dt::Float64 = 60.0,
    t_end::Float64 = 12 * 3600.0,
    Cs::Float64 = 2.0e5,
    Ks::Float64 = 2.0,
    T_deep::Float64 = 285.0,
    Rn0::Float64 = -50.0,
    rho_cp::Float64 = 1200.0,
    Ri_c::Float64 = 0.20,
    alpha::Float64 = 2.0,
    epsilon_c::Float64 = 0.05,
    eps_hyper::Float64 = 1e-3,
    l_m::Float64 = 15.0,
    S_eff::Float64 = 0.015
)
    return SEBModelConfig{Float64}(
        Nz, z_max, dt, t_end, Cs, Ks, T_deep, Rn0, rho_cp,
        Ri_c, alpha, epsilon_c, eps_hyper, l_m, S_eff
    )
end

@inline Phi_eps(x::T, eps_h::T) where {T<:AbstractFloat} = sqrt(x^2 + eps_h^2)

"""
    run_seb_simulation(use_cz0hr::Bool, config::SEBModelConfig)

Executes a 12-hour coupled SEB + 1D vertical diffusion integration.
"""
function run_seb_simulation(use_cz0hr::Bool, config::SEBModelConfig{T}) where {T<:AbstractFloat}
    # Grid Setup
    z = collect(range(1.5, config.z_max, length=config.Nz))
    dz = z - z
    g_grav = T(9.81)
    theta_0 = T(285.0)

    # State Vector Initialization: Initial uniform profile at 285.0 K
    theta = fill(T(285.0), config.Nz)
    Ts = T(285.0)

    Nt = round(Int, config.t_end / config.dt)
    time_hrs = zeros(T, Nt)
    Ts_history = zeros(T, Nt)
    H0_history = zeros(T, Nt)
    Kh_level1 = zeros(T, Nt)

    # Pre-allocate Tridiagonal Linear System
    A_tri = zeros(T, config.Nz, config.Nz)
    b_rhs = zeros(T, config.Nz)

    for n in 1:Nt
        t_sec = n * config.dt
        time_hrs[n] = t_sec / 3600.0

        # 1. Compute Raw Richardson Profile
        dtheta_dz = zeros(T, config.Nz)
        dtheta_dz = (theta - Ts) / z
        for i in 2:config.Nz
            dtheta_dz[i] = (theta[i] - theta[i-1]) / dz
        end

        Ri_raw = [(g_grav / theta_0) * dt_i / (config.S_eff^2) for dt_i in dtheta_dz]

        # 2. Evaluate Diffusivity Kh
        Kh = zeros(T, config.Nz)
        if !use_cz0hr
            # Unregularized C^0 Closure
            for i in 1:config.Nz
                g_val = one(T) - Ri_raw[i] / config.Ri_c
                Sh = (max(zero(T), g_val))^2
                Kh[i] = (config.l_m^2) * config.S_eff * Sh
            end
        else
            # C-Z0HR C^∞ Soft-Plus Closure
            d2Ri_dz2 = zeros(T, config.Nz)
            for i in 2:(config.Nz-1)
                d2Ri_dz2[i] = (Ri_raw[i+1] - 2*Ri_raw[i] + Ri_raw[i-1]) / (dz^2)
            end
            d2Ri_dz2 = d2Ri_dz2
            d2Ri_dz2[end] = d2Ri_dz2[end-1]

            for i in 1:config.Nz
                x_val = Ri_raw[i] - config.Ri_c
                Phi = Phi_eps(x_val, config.eps_hyper)
                C_scale = T(0.5) * abs(d2Ri_dz2[i]) * (dz^2) + config.epsilon_c * config.Ri_c

                Ri_reg = config.Ri_c + x_val * ((Phi + C_scale) / (Phi + (one(T) + config.alpha) * C_scale))

                g_reg = one(T) - Ri_reg / config.Ri_c
                softplus_val = T(0.5) * (g_reg + Phi_eps(g_reg, config.eps_hyper))
                Sh_reg = softplus_val^2
                Kh[i] = (config.l_m^2) * config.S_eff * Sh_reg
            end
        end

        # 3. SEB Backward-Euler Update for Ts
        gamma_coupling = config.rho_cp * Kh / z
        beta_denom = (config.Cs / config.dt) + gamma_coupling + config.Ks

        Ts_next = ( (config.Cs / config.dt) * Ts + config.Rn0 + gamma_coupling * theta + config.Ks * config.T_deep ) / beta_denom
        Ts = Ts_next

        # Downward Sensible Heat Flux H0 [W/m²]
        H0 = -config.rho_cp * Kh * (theta - Ts) / z

        # 4. Implicit 1D Vertical Column Diffusion for Potential Temperature
        A_tri .= zero(T)

        # Boundary level 1 (coupled to Ts)
        r_1 = Kh * config.dt / (dz * z)
        A_tri = one(T) + r_1
        A_tri = -r_1
        b_rhs = theta + (config.dt / (config.rho_cp * z)) * (config.rho_cp * Kh * Ts / z)

        # Interior levels
        for i in 2:(config.Nz-1)
            r_i = Kh[i] * config.dt / (dz^2)
            A_tri[i, i-1] = -r_i
            A_tri[i, i]   = one(T) + 2*r_i
            A_tri[i, i+1] = -r_i
            b_rhs[i]      = theta[i]
        end

        # Top boundary (no-flux)
        r_N = Kh[end] * config.dt / (dz^2)
        A_tri[end, end-1] = -r_N
        A_tri[end, end]   = one(T) + r_N
        b_rhs[end]        = theta[end]

        theta = A_tri \ b_rhs

        # Log Metrics
        Ts_history[n] = Ts
        H0_history[n] = H0
        Kh_level1[n]  = Kh
    end

    return time_hrs, Ts_history, H0_history, Kh_level1
end

function run_seb_benchmark()
    config = SEBModelConfig()
    println("Running 12-Hour SEB Cooling Simulations...")

    t_hrs, Ts_unreg, H0_unreg, Kh_unreg = run_seb_simulation(false, config)
    _,     Ts_cz0hr, H0_cz0hr, Kh_cz0hr = run_seb_simulation(true, config)

    @printf("\n=========================================================================\n")
    @printf("          SURFACE ENERGY BUDGET (SEB) COOLING CYCLE AUDIT (12 HOURS)\n")
    @printf("=========================================================================\n")
    @printf("Scheme         | Final Ts (K) | Cooling ΔTs (K) | Final H0 (W/m²) | Final Kh (m²/s)\n")
    @printf("-------------------------------------------------------------------------\n")
    @printf("Unregularized  | %12.2f | %15.2f | %15.2f | %15.2e\n",
            Ts_unreg[end], Ts_unreg[end] - 285.0, H0_unreg[end], Kh_unreg[end])
    @printf("C-Z0HR Reg.    | %12.2f | %15.2f | %15.2f | %15.2e\n",
            Ts_cz0hr[end], Ts_cz0hr[end] - 285.0, H0_cz0hr[end], Kh_cz0hr[end])
    @printf("=========================================================================\n")
end

run_seb_benchmark()