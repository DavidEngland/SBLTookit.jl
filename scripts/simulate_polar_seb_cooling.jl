#!/usr/bin/env julia
# =============================================================================
# ENHANCED POLAR 1D SEB SIMULATION WITH SOIL MOISTURE FREEZING & LATENT HEAT
# =============================================================================

using Printf
using LinearAlgebra

struct PolarSEBConfig{T<:AbstractFloat}
    # Atmospheric Grid & Time Parameters
    Nz::Int           # Number of vertical grid levels (default: 50)
    z_max::T          # Top of column [m] (default: 200.0)
    dt::T             # Timestep [s] (default: 60.0)
    t_end::T          # Total integration time [s] (default: 12 * 3600.0)
    
    # Surface Energy Budget Parameters
    Cs_base::T        # Base dry soil heat capacity [J/m²/K] (1.5e5)
    Ks::T             # Deep soil conductivity coefficient [W/m²/K] (1.5)
    T_deep::T         # Deep soil temperature [K] (275.0)
    Rn0::T            # Net surface radiative cooling [W/m²] (-45.0)
    rho_cp::T         # Air heat capacity constant [J/m³/K] (1280.0)
    rho::T            # Air density [kg/m³] (1.275)
    Ls::T             # Latent heat of sublimation/condensation [J/kg] (2.83e6)
    Lf::T             # Latent heat of fusion (freezing) [J/kg] (3.33e5)
    
    # Soil Moisture Freezing Parameters
    T_m::T            # Freezing point [K] (273.15)
    sigma_T::T        # Phase-change transition window width [K] (0.5)
    w_water::T        # Soil liquid water content [kg/m²] (4.0)
    
    # Closure Parameters
    Ri_c::T           # Critical Richardson number (0.20)
    alpha::T          # C-Z0HR damping control parameter (2.0)
    epsilon_c::T      # Curvature softening cap (0.05)
    eps_hyper::T      # C^∞ hyperbolic smoothing width (1e-3)
    kappa::T          # von Kármán constant (0.40)
    l_inf::T          # Asymptotic mixing length [m] (15.0)
    S_eff::T          # Ambient wind shear [s⁻¹] (0.045)
end

function PolarSEBConfig(;
    Nz::Int = 50,
    z_max::Float64 = 200.0,
    dt::Float64 = 60.0,
    t_end::Float64 = 12 * 3600.0,
    Cs_base::Float64 = 1.5e5,
    Ks::Float64 = 1.5,
    T_deep::Float64 = 275.0,
    Rn0::Float64 = -45.0,
    rho_cp::Float64 = 1280.0,
    rho::Float64 = 1.275,
    Ls::Float64 = 2.83e6,
    Lf::Float64 = 3.33e5,
    T_m::Float64 = 273.15,
    sigma_T::Float64 = 0.5,
    w_water::Float64 = 4.0,
    Ri_c::Float64 = 0.20,
    alpha::Float64 = 2.0,
    epsilon_c::Float64 = 0.05,
    eps_hyper::Float64 = 1e-3,
    kappa::Float64 = 0.40,
    l_inf::Float64 = 15.0,
    S_eff::Float64 = 0.045
)
    return PolarSEBConfig{Float64}(
        Nz, z_max, dt, t_end, Cs_base, Ks, T_deep, Rn0, rho_cp, rho, Ls, Lf,
        T_m, sigma_T, w_water, Ri_c, alpha, epsilon_c, eps_hyper, kappa, l_inf, S_eff
    )
end

@inline Phi_eps(x::T, eps_h::T) where {T<:AbstractFloat} = sqrt(x^2 + eps_h^2)

# Tetens saturation specific humidity [kg/kg] over ice/water
function q_sat(T_k::T) where {T<:AbstractFloat}
    T_c = T_k - T(273.15)
    e_sat = T(611.2) * exp((T(17.67) * T_c) / (T_c + T(243.5)))
    p_surf = T(101325.0)
    return T(0.622) * e_sat / p_surf
end

"""
    run_polar_seb_simulation(use_cz0hr, enable_freezing, enable_latent, config)

Executes a 12-hour coupled polar SEB + 1D vertical diffusion integration with
soil moisture phase-change latent heat release and frost deposition fluxes.
"""
function run_polar_seb_simulation(
    use_cz0hr::Bool,
    enable_freezing::Bool,
    enable_latent::Bool,
    config::PolarSEBConfig{T}
) where {T<:AbstractFloat}

    # Grid Setup
    z = collect(range(T(1.5), config.z_max, length=config.Nz))
    dz = z[2] - z[1]
    g_grav = T(9.81)
    theta_0 = T(275.0)
    
    # Height-dependent mixing length l_m(z)
    l_m = [(config.kappa * zi) / (one(T) + (config.kappa * zi) / config.l_inf) for zi in z]
    
    # State Vectors
    theta = fill(T(275.0), config.Nz)  # Initial air temperature [K]
    q_air = fill(T(0.0038), config.Nz) # Initial specific humidity [kg/kg]
    Ts = T(275.0)                      # Initial surface temperature [K]
    
    Nt = round(Int, config.t_end / config.dt)
    time_hrs = zeros(T, Nt)
    Ts_history = zeros(T, Nt)
    H0_history = zeros(T, Nt)
    LE0_history = zeros(T, Nt)
    Gs_history = zeros(T, Nt)
    Cs_eff_history = zeros(T, Nt)
    Kh_level1 = zeros(T, Nt)
    
    A_tri = zeros(T, config.Nz, config.Nz)
    b_rhs = zeros(T, config.Nz)
    
    for n in 1:Nt
        t_sec = n * config.dt
        time_hrs[n] = t_sec / T(3600.0)
        
        # 1. Temperature Gradients
        dtheta_dz = zeros(T, config.Nz)
        dtheta_dz[1] = (theta[1] - Ts) / z[1]
        for i in 2:config.Nz
            dtheta_dz[i] = (theta[i] - theta[i-1]) / dz
        end
        
        Ri_raw = [(g_grav / theta_0) * dtheta_dz[i] / (config.S_eff^2) for i in 1:config.Nz]
        
        # 2. Evaluate Diffusivity Kh
        Kh = zeros(T, config.Nz)
        if !use_cz0hr
            # Unregularized C^0 Closure
            for i in 1:config.Nz
                g_val = one(T) - Ri_raw[i] / config.Ri_c
                Sh = (max(zero(T), g_val))^2
                Kh[i] = (l_m[i]^2) * config.S_eff * Sh
            end
        else
            # C-Z0HR C^∞ Soft-Plus Closure
            d2Ri_dz2 = zeros(T, config.Nz)
            for i in 2:(config.Nz-1)
                d2Ri_dz2[i] = (Ri_raw[i+1] - 2*Ri_raw[i] + Ri_raw[i-1]) / (dz^2)
            end
            d2Ri_dz2[1] = d2Ri_dz2[2]
            d2Ri_dz2[end] = d2Ri_dz2[end-1]
            
            for i in 1:config.Nz
                x_val = Ri_raw[i] - config.Ri_c
                Phi = Phi_eps(x_val, config.eps_hyper)
                C_scale = T(0.5) * abs(d2Ri_dz2[i]) * (dz^2) + config.epsilon_c * config.Ri_c
                
                Ri_reg = config.Ri_c + x_val * ((Phi + C_scale) / (Phi + (one(T) + config.alpha) * C_scale))
                
                g_reg = one(T) - Ri_reg / config.Ri_c
                softplus_val = T(0.5) * (g_reg + Phi_eps(g_reg, config.eps_hyper))
                Sh_reg = softplus_val^2
                Kh[i] = (l_m[i]^2) * config.S_eff * Sh_reg
            end
        end
        
        Kh1 = Kh[1]
        
        # 3. Apparent Soil Heat Capacity C_s(T_s) with Soil Freezing Latent Heat Release
        if enable_freezing
            latent_freeze_cap = (config.w_water * config.Lf / (sqrt(T(2*pi)) * config.sigma_T)) *
                                exp(-T(0.5) * ((Ts - config.T_m) / config.sigma_T)^2)
            Cs_eff = config.Cs_base + latent_freeze_cap
        else:
            Cs_eff = config.Cs_base
        end
        Cs_eff_history[n] = Cs_eff
        
        # 4. Latent Heat Flux LE0 Coupling
        qs_surf = q_sat(Ts)
        gamma_q = enable_latent ? (config.rho * config.Ls * Kh1 / z[1]) : zero(T)
        gamma_h = config.rho_cp * Kh1 / z[1]
        
        # Linearize q_sat(Ts)
        T_c = Ts - T(273.15)
        dqs_dT = qs_surf * (T(17.67) * T(243.5)) / ((T_c + T(243.5))^2)
        
        # Backward-Euler SEB Equation
        A_seb = (Cs_eff / config.dt) + gamma_h + gamma_q * dqs_dT + config.Ks
        B_seb = (Cs_eff / config.dt) * Ts + config.Rn0 + gamma_h * theta[1] + gamma_q * (q_air[1] - qs_surf + dqs_dT * Ts) + config.Ks * config.T_deep
        
        Ts_next = B_seb / A_seb
        Ts = Ts_next
        
        # Re-evaluate Fluxes
        H0 = -config.rho_cp * Kh1 * (theta[1] - Ts) / z[1]
        LE0 = -config.rho * config.Ls * Kh1 * (q_air[1] - q_sat(Ts)) / z[1]
        Gs = config.Ks * (Ts - config.T_deep)
        
        # 5. Implicit Column Diffusion for theta
        A_tri .= zero(T)
        r1 = Kh1 * config.dt / (dz * z[1])
        A_tri[1, 1] = one(T) + r1
        A_tri[1, 2] = -r1
        b_rhs[1] = theta[1] + (config.dt / (config.rho_cp * z[1])) * (config.rho_cp * Kh1 * Ts / z[1])
        
        for i in 2:(config.Nz-1)
            ri = Kh[i] * config.dt / (dz^2)
            A_tri[i, i-1] = -ri
            A_tri[i, i]   = one(T) + 2*ri
            A_tri[i, i+1] = -ri
            b_rhs[i]      = theta[i]
        end
        
        rN = Kh[end] * config.dt / (dz^2)
        A_tri[end, end-1] = -rN
        A_tri[end, end]   = one(T) + rN
        b_rhs[end]        = theta[end]
        
        theta = A_tri \ b_rhs
        
        # Log History
        Ts_history[n] = Ts
        H0_history[n] = H0
        LE0_history[n] = LE0
        Gs_history[n] = Gs
        Kh_level1[n]  = Kh1
    end
    
    return time_hrs, Ts_history, H0_history, LE0_history, Gs_history, Cs_eff_history, Kh_level1
end

function run_polar_seb_benchmark()
    config = PolarSEBConfig()
    println("Running Humid Polar SEB Cooling Simulations (12 Hours)...")
    
    t_hrs, Ts_unreg, H0_unreg, LE0_unreg, Gs_unreg, Cs_unreg, Kh_unreg = run_polar_seb_simulation(false, true, true, config)
    _,     Ts_reg,   H0_reg,   LE0_reg,   Gs_reg,   Cs_reg,   Kh_reg   = run_polar_seb_simulation(true, true, true, config)
    _,     Ts_nofrz, H0_nofrz, LE0_nofrz, Gs_nofrz, Cs_nofrz, Kh_nofrz = run_polar_seb_simulation(true, false, true, config)
    
    @printf("\n=========================================================================================================\n")
    @printf("            HUMID POLAR SEB COOLING BENCHMARK: SOIL FREEZING & LATENT HEAT RELEASE\n")
    @printf("=========================================================================================================\n")
    @printf("Simulation Variant           | Final Ts (°C) | Final H0 (W/m²) | Final LE0 (W/m²) | Final Kh (m²/s)\n")
    @printf("---------------------------------------------------------------------------------------------------------\n")
    @printf("Unregularized (Quenched)      | %13.2f | %15.2f | %16.2f | %15.2e\n",
            Ts_unreg[end] - 273.15, H0_unreg[end], LE0_unreg[end], Kh_unreg[end])
    @printf("C-Z0HR (No Soil Freezing)     | %13.2f | %15.2f | %16.2f | %15.2e\n",
            Ts_nofrz[end] - 273.15, H0_nofrz[end], LE0_nofrz[end], Kh_nofrz[end])
    @printf("C-Z0HR Full Polar SEB        | %13.2f | %15.2f | %16.2f | %15.2e\n",
            Ts_reg[end] - 273.15, H0_reg[end], LE0_reg[end], Kh_reg[end])
    @printf("=========================================================================================================\n\n")
end

run_polar_seb_benchmark()
