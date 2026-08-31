using Test
using NCDatasets
using Dates
using LinearAlgebra
using Statistics

if !isdefined(Main, :NetCDFIngestionEngine)
    include(joinpath(@__DIR__, "..", "src", "netcdf-ingestion-engine-v2.jl"))
end
using .NetCDFIngestionEngine

@testset "NetCDFIngestionEngine.jl Suite" begin

    @testset "Helper Utilities & Sentinels" begin
        @test NetCDFIngestionEngine.nanmax_abs([1.0, -5.5, NaN, 3.2]) == 5.5
        @test NetCDFIngestionEngine.nanmax_abs([NaN, missing]) == -Inf

        @test NetCDFIngestionEngine.looks_like_celsius([-10.0, 0.0, 25.0]) == true
        @test NetCDFIngestionEngine.looks_like_celsius([263.15, 273.15, 298.15]) == false

        @test NetCDFIngestionEngine.looks_like_julian_day([1.5, 180.2, 365.0]) == true
        @test NetCDFIngestionEngine.looks_like_julian_day([1000.0, 2000.0]) == false

        raw_data = Any[12.5, 999.0, -9999.0, missing, NaN, 99999.0]
        scrubbed = NetCDFIngestionEngine.scrub_sentinels!(copy(raw_data))
        @test scrubbed[1] == 12.5
        @test isnan(scrubbed[2])
        @test isnan(scrubbed[3])
        @test isnan(scrubbed[4])
        @test isnan(scrubbed[5])
        @test isnan(scrubbed[6])

        vec_t = [1.0, 2.0, 3.0]
        mat_replicated = NetCDFIngestionEngine.ensure_2d_matrix(vec_t, 2, 3)
        @test size(mat_replicated) == (2, 3)
        @test mat_replicated[1, :] == vec_t
        @test mat_replicated[2, :] == vec_t
    end

    @testset "Tower Suffix Ingestion (Mode A)" begin
        mktempdir() do tmpdir
            nc_path = joinpath(tmpdir, "test_tower.nc")

            Dataset(nc_path, "c") do ds
                defDim(ds, "time", 2)

                v_time = defVar(ds, "time", Float64, ("time",))
                v_time[:] = [1.0, 2.0]

                v_u10 = defVar(ds, "u_10m", Float64, ("time",))
                v_v10 = defVar(ds, "v_10m", Float64, ("time",))
                v_t10 = defVar(ds, "tc_10m", Float64, ("time",))

                v_u20 = defVar(ds, "u_20m", Float64, ("time",))
                v_v20 = defVar(ds, "v_20m", Float64, ("time",))
                v_t20 = defVar(ds, "tc_20m", Float64, ("time",))

                v_u10[:] = [3.0, 4.0]
                v_v10[:] = [4.0, 3.0]
                v_t10[:] = [15.0, 16.0]

                v_u20[:] = [6.0, 8.0]
                v_v20[:] = [8.0, 6.0]
                v_t20[:] = [10.0, 11.0]
            end

            campaign = ingest_netcdf_gspt(nc_path)

            @test campaign isa IngestedCampaign
            @test length(campaign.profiles) == 2

            prof1 = campaign.profiles[1]
            @test prof1.z == [10.0, 20.0]
            @test prof1.U ≈ [5.0, 10.0]
            @test prof1.theta[1] ≈ (15.0 + 273.15 + 0.0098 * 10.0) atol = 1e-4
        end
    end

    @testset "Explicit Profile Ingestion & Grid Sorting (Mode B)" begin
        mktempdir() do tmpdir
            nc_path = joinpath(tmpdir, "test_profile.nc")

            Dataset(nc_path, "c") do ds
                defDim(ds, "z", 3)
                defDim(ds, "time", 2)

                v_z = defVar(ds, "z", Float64, ("z",))
                v_z[:] = [50.0, 20.0, 2.0]

                v_time = defVar(ds, "time", Float64, ("time",))
                v_time[:] = [100.0, 101.0]

                v_u = defVar(ds, "u", Float64, ("z", "time"))
                v_v = defVar(ds, "v", Float64, ("z", "time"))
                v_t = defVar(ds, "theta", Float64, ("z", "time"))

                v_u[:, :] = [5.0 5.0; 3.0 3.0; 1.0 1.0]
                v_v[:, :] = [0.0 0.0; 0.0 0.0; 0.0 0.0]
                v_t[:, :] = [300.0 301.0; 295.0 296.0; 290.0 291.0]
            end

            campaign = ingest_netcdf_gspt(nc_path)
            prof = campaign.profiles[1]

            @test prof.z == [2.0, 20.0, 50.0]
            @test prof.U == [1.0, 3.0, 5.0]
            @test prof.theta == [290.0, 295.0, 300.0]
        end
    end

    @testset "GABLS3 Benchmark Ingestion" begin
        mktempdir() do tmpdir
            nc_path = joinpath(tmpdir, "gabls3_scm.nc")

            Dataset(nc_path, "c") do ds
                defDim(ds, "zf", 3)
                defDim(ds, "time", 2)

                v_zf = defVar(ds, "zf", Float64, ("zf",))
                v_zf[:] = [10.0, 50.0, 100.0]

                v_time = defVar(ds, "time", Float64, ("time",))
                v_time[:] = [0.0, 1.0]

                v_u = defVar(ds, "u", Float64, ("zf", "time"))
                v_v = defVar(ds, "v", Float64, ("zf", "time"))
                v_t = defVar(ds, "ta", Float64, ("zf", "time"))

                v_u[:, :] = [2.0 3.0; 4.0 5.0; 6.0 7.0]
                v_v[:, :] = [0.0 0.0; 0.0 0.0; 0.0 0.0]
                v_t[:, :] = [10.0 11.0; 12.0 13.0; 14.0 15.0]
            end

            scm_data = ingest_gabls3_netcdf(nc_path)

            @test scm_data.Nz == 3
            @test scm_data.Nt == 2
            @test scm_data.z == [10.0, 50.0, 100.0]
            @test scm_data.theta[1, 1] ≈ (10.0 + 273.15 + 0.0098 * 10.0) atol = 1e-4
            @test scm_data.U[1, 1] ≈ 2.0
        end
    end
end
