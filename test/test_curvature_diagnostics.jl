#!/usr/bin/env julia
# SBLToolKit.jl: Curvature Diagnostics Test Suite
# test/test_curvature_diagnostics.jl`
using Test
using Statistics
using Random

# Load local standalone module without modifying package exports
if !@isdefined(CurvatureDiagnostics)
    include(joinpath(@__DIR__, "..", "src", "CurvatureDiagnostics.jl"))
end
using .CurvatureDiagnostics

@testset "CurvatureDiagnostics Diagnostic Suite" begin

    @testset "1. Output Struct & Vector Bounding" begin
        Nz = 100
        z = collect(range(1.0, stop=100.0, length=Nz))
        u_raw = 0.1 .* z
        v_raw = zeros(Nz)
        theta_raw = 300.0 .+ 0.05 .* z

        out = process_dns_profile(z, u_raw, v_raw, theta_raw)

        @test out isa CurvatureOutput
        @test length(out.z) == Nz
        @test length(out.Rig) == Nz
        @test length(out.zeta) == Nz
        @test length(out.zeta_z) == Nz
        @test length(out.zeta_zz) == Nz
        @test length(out.C_M) == Nz
        @test length(out.classification) == Nz

        # Assert bounded mapping fraction [0, 1]
        @test all((0.0 .<= out.C_M) .& (out.C_M .<= 1.0))
    end

    @testset "2. Linear Baseline Profile (No Coordinate Fold)" begin
        Nz = 100
        z = collect(range(1.0, stop=100.0, length=Nz))
        u_raw = 0.2 .* z
        v_raw = zeros(Nz)
        theta_raw = 300.0 .+ 0.02 .* z

        out = process_dns_profile(z, u_raw, v_raw, theta_raw)

        # Constant shear and linear stratification produce no coordinate fold
        @test !any(out.classification .== :PureCoordinateFold)
    end

    @testset "3. Synthetic Jet Fold Detection & Spatial Localization" begin
        Nz = 150
        z = collect(range(2.0, stop=300.0, length=Nz))

        # Synthetic LLJ profile with zero shear at z = 150 m
        u_raw = @. 10.0 * sin(π * z / 300.0)
        v_raw = zeros(Nz)
        theta_raw = @. 300.0 + 5.0 * (z / 300.0)^2

        out = process_dns_profile(z, u_raw, v_raw, theta_raw)

        # Enforce enum symbol validity
        valid_symbols = Set([:PureCoordinateFold, :PureDynamicFold, :Ambiguous, :HybridResonantFold])
        @test all(s -> s in valid_symbols, out.classification)

        # Enforce active detection of the zero-shear coordinate fold
        @test any(out.classification .== :PureCoordinateFold)

        # Verify spatial localization within 15 m of jet nose (z = 150 m)
        fold_indices = findall(s -> s == :PureCoordinateFold, out.classification)
        fold_heights = z[fold_indices]
        @test any(abs.(fold_heights .- 150.0) .<= 15.0)
    end

    @testset "4. Noise Robustness & GCV Smoothing" begin
        Random.seed!(42)
        Nz = 100
        z = collect(range(1.0, stop=100.0, length=Nz))
        u_clean = 0.2 .* z
        v_clean = zeros(Nz)
        theta_clean = 300.0 .+ 0.02 .* z

        # Gaussian perturbation
        u_noisy = u_clean .+ 0.01 .* randn(Nz)
        v_noisy = v_clean .+ 0.01 .* randn(Nz)
        theta_noisy = theta_clean .+ 0.01 .* randn(Nz)

        out_clean = process_dns_profile(z, u_clean, v_clean, theta_clean)
        out_noisy = process_dns_profile(z, u_noisy, v_noisy, theta_noisy)

        # GCV spline filter suppresses derivative noise blowup in C_M
        @test mean(abs.(out_noisy.C_M .- out_clean.C_M)) < 0.15
    end

end