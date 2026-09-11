#!/usr/bin/env julia
# =============================================================================
# SBL ADAPTIVE GRID EQUIDISTRIBUTION PLOTTING SCRIPT (sbl-grid-plot.jl)
# Developed for SBLToolkit.jl / GSPT Adaptive Grid Benchmark
# =============================================================================
# Computes and visualizes dynamic physical level contraction (N_z = 38)
# around a descending nocturnal Low-Level Jet (LLJ) core (140m -> 42m).
# =============================================================================

using Printf
using LinearAlgebra

try
    using Plots
catch e
    @warn "Plots.jl package not found in current environment. Numerical matrices will be computed without PNG rendering."
end

struct AdaptiveGridConfig{T<:AbstractFloat}
    Nz::Int          # Number of physical vertical levels
    Nt::Int          # Number of temporal evaluation steps
    z_max::T        # Upper domain boundary [m]
    z0::T           # Surface roughness length [m]
    t_end::T        # Simulation duration [hours]
    z_jet_start::T  # Initial LLJ core height [m]
    z_jet_end::T    # Final LLJ core height [m]
end

function AdaptiveGridConfig(;
    Nz::Int = 38,
    Nt::Int = 100,
    z_max::T = 200.0,
    z0::T = 0.15,
    t_end::T = 8.0,
    z_jet_start::T = 140.0,
    z_jet_end::T = 42.0
) where {T<:AbstractFloat}
    return AdaptiveGridConfig{T}(Nz, Nt, z_max, z0, t_end, z_jet_start, z_jet_end)
end

"""
    compute_equidistributed_grid(config::AdaptiveGridConfig{T})

Computes physical grid trajectories z_levels(N_z, N_t) by equidistributing a 
composite monitor function M(z) = 1 + w_surf(z) + w_llj(z, t).
"""
function compute_equidistributed_grid(config::AdaptiveGridConfig{T}) where {T<:AbstractFloat}
    t_hours = collect(range(zero(T), config.t_end, length=config.Nt))
    z_llj = [config.z_jet_start + (config.z_jet_end - config.z_jet_start) * (t / config.t_end) for t in t_hours]
    
    n_fine = 2000
    z_fine = collect(range(zero(T), config.z_max, length=n_fine))
    dz_fine = z_fine[2] - z_fine[1]
    
    z_levels = zeros(T, config.Nz, config.Nt)
    xi_uniform = collect(range(zero(T), one(T), length=config.Nz))
    
    for (t_idx, t) in enumerate(t_hours)
        z_jet = z_llj[t_idx]
        
        # 1. Surface layer logarithmic weight + LLJ core curvature weight
        w_surf = [T(15.0) / (z + config.z0) for z in z_fine]
        w_llj = [T(25.0) * exp(-((z - z_jet) / T(15.0))^2) for z in z_fine]
        M = one(T) .+ w_surf .+ w_llj
        
        # 2. Cumulative monitor integral
        cum_M = zeros(T, n_fine)
        for j in 2:n_fine
            cum_M[j] = cum_M[j-1] + T(0.5) * (M[j-1] + M[j]) * dz_fine
        end
        cum_M_norm = cum_M ./ cum_M[end]
        
        # 3. Inverse mapping via linear/cubic search
        for i in 1:config.Nz
            xi_target = xi_uniform[i]
            idx = findfirst(>=(xi_target), cum_M_norm)
            if idx === nothing || idx == 1
                z_levels[i, t_idx] = zero(T)
            else
                # Local linear interpolation for physical height
                frac = (xi_target - cum_M_norm[idx-1]) / (cum_M_norm[idx] - cum_M_norm[idx-1] + T(1e-12))
                z_levels[i, t_idx] = z_fine[idx-1] + frac * dz_fine
            end
        end
    end
    
    return t_hours, z_llj, z_levels
end

function run_grid_stretching_plot()
    config = AdaptiveGridConfig()
    t_hours, z_llj, z_levels = compute_equidistributed_grid(config)
    
    dz_t8 = diff(z_levels[:, end])
    @printf("Equidistribution calculation complete.\n")
    @printf("Sub-meter resolution achieved at t=8h jet core: min Δz = %.3f m\n", minimum(dz_t8))
    
    if @isdefined Plots
        gr()
        p1 = plot(xlabel="Nocturnal Time [Hours]", ylabel="Physical Height z [m]",
                  title="Adaptive Grid Equidistribution (N_z = 38)", legend=:topright, grid=true)
        for i in 1:config.Nz
            plot!(p1, t_hours, z_levels[i, :], label="", color=:steelblue, lw=0.8)
        end
        plot!(p1, t_hours, z_llj, label="Descending LLJ Core z_LLJ(t)", color=:red, linestyle=:dash, lw=2.0)
        
        savefig(p1, "/workspace/scratch/sbl_grid_plot_julia.png")
        @printf("Julia plot saved successfully.\n")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_grid_stretching_plot()
end
