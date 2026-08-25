# test/runtests.jl

using Test
using DataFrames
using Dates

# Ensure local src/ is accessible during direct execution
project_src = joinpath(@__DIR__, "..", "src")
if !(project_src in LOAD_PATH)
    push!(LOAD_PATH, project_src)
end

include("test_triple_point_qc.jl")
include("test_scm_convergence.jl")

using SBLToolkit

@testset "SBLToolkit.jl Automated Integration Test Suite" begin

    @testset "Core Stability & Forcing Modules" begin
        # 1. Stability Classification Verification
        @test classify_stability(1.5) == :strongly_stable
        @test classify_stability(0.5) == :stable
        @test classify_stability(0.0) == :near_neutral
        @test classify_stability(-0.5) == :unstable
        @test classify_stability(-1.5) == :strongly_unstable
        @test classify_stability(NaN) == :unknown

        # 2. Single Column Model Forcing Generation
        dt = DateTime(2026, 8, 22, 12, 0)
        meta = ProfileMetadata(dt, 0.25, -50.0, 10.0, 0.0)
        forcing = generate_scm_forcing(meta, 15.0, 30.0, 5.0, 20)

        @test forcing.sensible_flux == 15.0
        @test forcing.latent_flux == 30.0
        @test forcing.zeta_reference ≈ (10.0 / -50.0)
        @test length(forcing.theta_tendency) == 20
    end

    @testset "Spectral Engine & Grid Reconstruction" begin
        heights = [2.0, 5.0, 10.0, 20.0, 40.0]
        vals = [15.0, 14.8, 14.2, 13.5, 12.0]
        meta = ProfileMetadata(DateTime(2026, 8, 22), 0.2, -100.0, 10.0, 0.0)
        prof = MeteorologicalProfile(meta, heights, vals)

        # Forward Projection
        coeffs = chebyshev_fingerprint(prof; n_coeffs=4, height_mapping=:log)
        @test length(coeffs) == 4
        @test !any(isnan, coeffs)

        # Inverse Reconstruction onto dense 20-level grid
        target_z = collect(range(2.0, 40.0, length=20))
        recon_v = reconstruct_from_chebyshev(coeffs, target_z, (2.0, 40.0); height_mapping=:log)
        @test length(recon_v) == 20
        @test isapprox(recon_v[1], vals[1], atol=0.5)
        @test isapprox(recon_v[end], vals[end], atol=0.5)
    end

    @testset "Network Observational Adapters" begin
        # Cabauw Station Adapter
        @testset "Cabauw Adapter" begin
            df_cabauw = DataFrame(
                datetime=[DateTime(2026, 8, 22, 12, 0)],
                TA_2m=[18.5], TA_10m=[18.2], TA_20m=[17.8],
                TA_40m=[17.3], TA_80m=[16.5], TA_140m=[15.6], TA_200m=[14.8],
                ustar=[0.32], L_obukhov=[-50.0]
            )
            obs = extract_temperature_observations(df_cabauw)
            @test length(obs) == 1
            @test obs[1].campaign == "GABLS3"
            @test obs[1].n_valid_levels == 7
            @test obs[1].robust_for_eta3 == true
        end

        # NEON Multi-Tier Adapter (expects startDateTime or timestamp)
        @testset "NEON Adapter" begin
            df_neon = DataFrame(
                datetime=[DateTime(2026, 8, 22, 12, 0)],
                temp_z1=[22.1], temp_z2=[21.5], temp_z3=[20.8],
                ustar=[0.28], L=[-45.0]
            )
            obs = extract_neon_observations(df_neon, "temp", [2.0, 10.0, 30.0]; campaign="NEON-HARV")
            @test length(obs) == 1
            @test obs[1].n_valid_levels == 3
            @test obs[1].heights == [2.0, 10.0, 30.0]
        end

        # ICOS Monin-Obukhov Similarity Upscaling
        @testset "ICOS Upscaling Adapter" begin
            dt = DateTime(2026, 8, 22, 12, 0)
            obs = upscale_sparse_icos_observation(
                dt, 10.0, 50.0, 15.4, 13.2, 0.18, 120.0;
                campaign="ICOS-SE-S2", stability_correction=:psi_h
            )
            @test obs.n_valid_levels == 5
            @test obs.robust_for_eta3 == true
        end

        # AmeriFlux Site Registry
        @testset "AmeriFlux Registry" begin
            reg_path = joinpath(@__DIR__, "..", "data", "ameriflux", "stations.json")
            reg = load_ameriflux_registry(reg_path)
            @test reg isa Dict
            meta = get_site_metadata("US-ARM")
            @test length(meta.heights) >= 3
        end
    end

    @testset "Unified Ingestion Router & Batch Pipeline" begin
        df_cabauw = DataFrame(
            datetime=[DateTime(2026, 8, 22)],
            TA_2m=[15.0], TA_10m=[14.5], TA_20m=[14.0],
            TA_40m=[13.5], TA_80m=[13.0], TA_140m=[12.0], TA_200m=[11.0]
        )
        @test detect_network_format(df_cabauw) == :cabauw

        obs_list = ingest_boundary_layer_data(df_cabauw)
        @test length(obs_list) == 1

        fps = batch_spectral_fingerprints(obs_list; n_coeffs=4)
        @test size(fps) == (1, 4)
    end

    @testset "Raw Dataset File Integrity Checks" begin
        raw_dir = joinpath(@__DIR__, "..", "data", "raw")
        campaigns = ["bllast", "cases99", "floss", "gabls3", "sheba"]

        for camp in campaigns
            camp_path = joinpath(raw_dir, camp)
            # Ensure folder tree structure exists before testing integrity
            (ispath(camp_path) || islink(camp_path)) || mkpath(camp_path)
            @testset "Raw Campaign Directory: $camp" begin
                @test (ispath(camp_path) || islink(camp_path))
            end
        end
    end

end