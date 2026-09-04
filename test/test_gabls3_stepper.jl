using Test
using SBLToolkit

@testset "GABLS3 implicit stepper" begin
    z = collect(5.0:5.0:30.0)
    config = SCMConfig(z; dt=60.0, theta_flux=-0.01, km_background=0.001, kh_background=0.001)
    gate_params = BifurcationGatingParams(epsilon_on=-0.10, epsilon_off=-0.05, zeta_z_tol=1.0)
    gs_config = GSPTModelConfig()

    baseline = initialize_cabauw_state(z)
    initial_theta = copy(baseline.theta_v)
    diagnostics = step_scm!(baseline, config, gate_params, gs_config; enable_gating=false)
    @test diagnostics.gated_levels == 0
    @test all(isfinite, baseline.theta_v)
    @test all(isfinite, baseline.u)
    @test baseline.theta_v != initial_theta
    @test isfinite(diagnostics.t2m)
    @test isfinite(diagnostics.boundary_layer_height)
    @test isfinite(diagnostics.surface_sensible_heat_flux)

    slab_config = SCMConfig(z; dt=60.0, surface_mode=:slab, radiation_max=400.0,
        longwave_loss=80.0, radiation_peak_time=43_200.0, soil_temperature=279.0,
        ground_conductance=2.0)
    slab_state = initialize_cabauw_state(z)
    initial_skin_temperature = slab_state.skin_temperature
    slab_diagnostics = step_scm!(slab_state, slab_config, gate_params, gs_config; enable_gating=false)
    @test net_radiation(slab_config, 43_200.0) == 320.0
    @test net_radiation(slab_config, 0.0) == -80.0
    @test slab_diagnostics.net_radiation == -80.0
    @test slab_diagnostics.surface_sensible_heat_flux < 0.0
    @test isfinite(slab_diagnostics.ground_heat_flux)
    @test slab_state.skin_temperature ≈ initial_skin_temperature + slab_config.dt *
        (slab_diagnostics.net_radiation - slab_diagnostics.surface_sensible_heat_flux -
         slab_diagnostics.ground_heat_flux) / slab_config.surface_heat_capacity
    @test slab_state.elapsed_seconds == slab_config.dt

    attracting_state = SCMState(z, fill(2.0, length(z)), zeros(length(z)),
        fill(280.0, length(z)), fill(0.01, length(z)), 279.0)
    gated_diagnostics = step_scm!(attracting_state, config, gate_params, gs_config; enable_gating=true)
    @test gated_diagnostics.gated_levels == length(z)
    @test all(gated_diagnostics.lambda_f .<= gate_params.epsilon_on)

    solution = zeros(3)
    solve_tridiagonal!(solution, [-1.0, -1.0], [2.0, 2.0, 2.0], [-1.0, -1.0], [1.0, 0.0, 1.0])
    @test solution ≈ [1.0, 1.0, 1.0]
end