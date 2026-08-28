using Test

@testset "Pluto notebook structure" begin
    notebook_path = joinpath(@__DIR__, "..", "scripts", "fast_slow_tke_pluto.jl")
    notebook = read(notebook_path, String)

    @test !contains(notebook, "00000000-0000-0000-0000-000000000001")
    @test !contains(notebook, "00000000-0000-0000-0000-000000000002")
    @test contains(notebook, "11111111-1111-1111-1111-111111111111")
    @test contains(notebook, "22222222-2222-2222-2222-222222222222")
end
