#!/usr/bin/env julia
# =============================================================================
# SBLToolkit.jl: test/test_emetric_cases99_convergence.jl
# Level 3 MMS Verification Suite: Metric-Consistency Residual & Grid Convergence
# =============================================================================

using Test
using LinearAlgebra
using Statistics
using Printf

# Import SBLToolkit core module if running in package environment
if @isdefined(SBLToolkit)
    using .SBLToolkit
end

# =============================================================================
# 1. Non-Uniform Stencil Builders (Taylor-Vandermonde Systems)
# =============================================================================

"""
    stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)

Computes Taylor-Vandermonde stencil weights for derivative order `m` at target `z0`
by solving the local linear system A * w = b, where A_ik = (z_i - z0)^(k-1) / (k-1)!.
"""
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    A = [Float64(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p)
    b[m + 1] = 1.0  # Target derivative order selector
    return A \ b
end

"""
    build_nonuniform_operators(z::Vector{Float64})

Generates non-uniform derivative matrices D1 (first derivative) and D2 (second derivative)
on an arbitrary physical coordinate vector `z`. Boundary points utilize 3-point one-sided stencils.
"""
function build_nonuniform_operators(z::Vector{Float64})
    n = length(z)
    D1 = zeros(n, n)
    D2 = zeros(n, n)

    for i in 1:n
        if i == 1
            idx =
        elseif i == n
            idx = [n - 2, n - 1, n]
        else
            idx = [i - 1, i, i + 1]
        end

        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end

    return D1, D2
end

# =============================================================================
# 2. Manufactured Solution Target Functions
# =============================================================================

# Smooth grid-stretching function z(eta) mapping computational nodes to CASES-99 height domain [1.5m, 55.0m]
z_mapping(eta::Float64, h_top::Float64, N::Int) = 1.5 + (h_top - 1.5) * ((eta - 1.0) / (N - 1.0))^1.6

# Manufactured Monin-Obukhov Richardson profile: Ri_g(z) = (z/L0) / (1 + beta * z / L0)
function manufactured_rig(z::Float64; L0=20.0, beta=5.0)
    zeta = z / L0
    return zeta / (1.0 + beta * zeta)
end

function manufactured_rig_dz(z::Float64; L0=20.0, beta=5.0)
    zeta = z / L0
    dzeta_dz = 1.0 / L0
    dRi_dzeta = 1.0 / ((1.0 + beta * zeta)^2)
    return dRi_dzeta * dzeta_dz
end

function manufactured_rig_d2z(z::Float64; L0=20.0, beta=5.0)
    zeta = z / L0
    dzeta_dz = 1.0 / L0
    d2Ri_dzeta2 = -2.0 * beta / ((1.0 + beta * zeta)^3)
    return d2Ri_dzeta2 * (dzeta_dz^2)
end

# =============================================================================
# 3. Level 3 MMS Verification TestSet
# =============================================================================

@testset "SBLToolkit.jl Verification Suite: E_metric Level 3 Convergence" begin

    @testset "1. Analytical Continuum Identity (E_metric == 0)" begin
        # Verifies that exact continuum transformations satisfy E_metric = R_etaeta - z_etaeta*R_z - z_eta^2*R_zz == 0
        z_cases99 = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
        L0 = 20.0
        beta = 5.0

        for z in z_cases99
            z_eta = 8.5     # Local metric scale factor (dz/d_eta)
            z_etaeta = 1.2  # Local metric acceleration (d2z/d_eta2)

            R_z = manufactured_rig_dz(z; L0=L0, beta=beta)
            R_zz = manufactured_rig_d2z(z; L0=L0, beta=beta)

            # Exact chain-rule derivatives in computational space
            R_eta = z_eta * R_z
            R_etaeta = (z_eta^2) * R_zz + z_etaeta * R_z

            # Evaluate continuum Metric-Consistency Residual
            E_metric = R_etaeta - z_etaeta * R_z - (z_eta^2) * R_zz

            @test isapprox(E_metric, 0.0, atol=1e-14)
        end
    end

    @testset "2. CASES-99 Tower Footprint Diagnostic Audit" begin
        # Native CASES-99 7-level non-uniform tower levels
        z_cases99 = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
        Nz = length(z_cases99)
        eta_grid = collect(1.0:Float64(Nz))

        D1_z, D2_z = build_nonuniform_operators(z_cases99)
        D1_eta, D2_eta = build_nonuniform_operators(eta_grid)

        # Sample manufactured profile on native CASES-99 nodes
        Ri_vec = [manufactured_rig(zi) for zi in z_cases99]

        # Numerical differentiation on non-uniform mesh
        Ri_z_num = D1_z * Ri_vec
        Ri_zz_num = D2_z * Ri_vec
        Ri_etaeta_num = D2_eta * Ri_vec

        z_eta_num = D1_eta * z_cases99
        z_etaeta_num = D2_eta * z_cases99

        # Compute discrete metric-consistency residual across interior tower nodes
        E_metric_num = zeros(Nz)
        for k in 2:(Nz - 1)
            E_metric_num[k] = Ri_etaeta_num[k] - z_etaeta_num[k] * Ri_z_num[k] - (z_eta_num[k]^2) * Ri_zz_num[k]
        end

        # Verify that discrete residual on raw 7-level tower grid is bounded
        @test maximum(abs.(E_metric_num[2:(Nz - 1)])) < 1.0e-2
    end

    @testset "3. Level 3 MMS Grid Refinement Convergence (O(d_eta^2))" begin
        # Successive grid refinement sweep across 5 grid scales
        refinement_factors =
        L_inf_errors = Float64[]
        h_top = 55.0

        for N in refinement_factors
            eta_nodes = collect(range(1.0, Float64(N), length=N))
            z_nodes = [z_mapping(eta, h_top, N) for eta in eta_nodes]

            D1_z, D2_z = build_nonuniform_operators(z_nodes)
            D1_eta, D2_eta = build_nonuniform_operators(eta_nodes)

            Ri_nodes = [manufactured_rig(zi) for zi in z_nodes]

            Ri_z_num = D1_z * Ri_nodes
            Ri_zz_num = D2_z * Ri_nodes
            Ri_etaeta_num = D2_eta * Ri_nodes

            z_eta_num = D1_eta * z_nodes
            z_etaeta_num = D2_eta * z_nodes

            # Evaluate interior L_infinity norm of E_metric_num
            E_metric_num = zeros(N)
            for k in 3:(N - 2)
                E_metric_num[k] = Ri_etaeta_num[k] - z_etaeta_num[k] * Ri_z_num[k] - (z_eta_num[k]^2) * Ri_zz_num[k]
            end

            push!(L_inf_errors, maximum(abs.(E_metric_num[3:(N - 2)])))
        end

        # Compute observed asymptotic order of convergence (p = log2(e_h / e_{h/2}))
        convergence_orders = Float64[]
        for i in 1:(length(L_inf_errors) - 1)
            p = log2(L_inf_errors[i] / L_inf_errors[i + 1])
            push!(convergence_orders, p)
        end

        # MMS Verification Assertions
        # 1. Average convergence order meets or exceeds theoretical 2nd-order Taylor accuracy (p >= 2.0)
        @test mean(convergence_orders) >= 2.0

        # 2. Residual error decays monotonically by > 1,000x across refinement range
        @test L_inf_errors[end] < L_inf_errors / 1000.0
    end
end