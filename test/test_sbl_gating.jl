using Test
using SBLToolkit

@testset "SBL geometry-aware gating diagnostics" begin
    params = BifurcationGatingParams(
        epsilon_on=-0.10,
        epsilon_off=-0.05,
        zeta_z_tol=1e-2,
        km_floor=0.1,
        kh_floor=0.01,
    )

    @testset "Coordinate fold with attracting fast mode" begin
        state = GatingState(Float64)
        lambda_f = extract_fast_eigenvalue([-0.412 0.0; 0.0 -0.01], state)
        K_m, K_h, is_gated = update_gating_state!(state, lambda_f, 1e-4, 0.005, 0.001, params)

        @test lambda_f == -0.412
        @test is_gated
        @test K_m == params.km_floor
        @test K_h == params.kh_floor
    end

    @testset "Weak fast mode preserves diagnostic diffusivity" begin
        state = GatingState(Float64)
        state.is_active = true
        K_m, K_h, is_gated = update_gating_state!(state, -0.02, 1e-4, 0.005, 0.001, params)

        @test !state.is_active
        @test !is_gated
        @test K_m == 0.005
        @test K_h == 0.001
    end

    @testset "Hysteresis and independent geometry criterion" begin
        state = GatingState(Float64)
        update_gating_state!(state, -0.12, 0.0, 0.01, 0.001, params)
        @test state.is_active

        _, _, is_gated = update_gating_state!(state, -0.075, 0.0, 0.01, 0.001, params)
        @test is_gated

        K_m, K_h, is_gated = update_gating_state!(state, -0.075, 0.02, 0.01, 0.001, params)
        @test state.is_active
        @test !is_gated
        @test K_m == 0.01
        @test K_h == 0.001

        _, _, is_gated = update_gating_state!(state, -0.02, 0.0, 0.01, 0.001, params)
        @test !is_gated
        @test !state.is_active
    end

    @testset "Eigenvector continuation across ordering swaps" begin
        state = GatingState(Float64)
        lambda_1 = extract_fast_eigenvalue([-0.412 0.0; 0.0 -0.01], state)
        lambda_2 = extract_fast_eigenvalue([-0.01 0.0; 0.0 -0.412], state)

        @test lambda_1 == -0.412
        @test lambda_2 == -0.01
        @test_throws ArgumentError extract_fast_eigenvalue([0.0 -1.0; 1.0 0.0], state)
        @test extract_fast_eigenvalue([0.0 -1.0; 1.0 0.0], state; complex_policy=:real_part) == 0.0
        @test isnothing(state.prev_v_f)
    end

    @testset "Gate diagnostics audit the geometric override" begin
        continuity_params = BifurcationGatingParams(
            epsilon_on=-0.10, epsilon_off=-0.05, zeta_z_tol=1e-2,
            min_mode_overlap=0.8,
        )
        state = GatingState(Float64)
        diagnostics = evaluate_gate_step!(
            state, [-0.412 0.0; 0.0 -0.01], 1e-4, true, continuity_params,
        )
        @test diagnostics.lambda_f == -0.412
        @test diagnostics.mode_overlap == 1.0
        @test diagnostics.gate_active
        @test diagnostics.is_coordinate_regular
        @test diagnostics.q_override

        diagnostics = evaluate_gate_step!(
            state, [-0.211 -0.201; -0.201 -0.211], 1e-4, true, continuity_params,
        )
        @test diagnostics.mode_overlap ≈ inv(sqrt(2))
        @test !diagnostics.gate_active
        @test !diagnostics.q_override
    end

    @testset "GABLS3 Jacobian adapter and fold audit" begin
        config = GSPTModelConfig()
        jacobian = map_gabls3_to_jacobian(2, [0.1, 0.05], [0.2, 0.1], [0.01, 0.02], config)
        fold_ratio = compute_fold_ratio(1e-6, 0.1, 1e-4)

        @test size(jacobian) == (2, 2)
        @test all(isfinite, jacobian)
        @test 0.0 <= fold_ratio <= 1.0
        @test fold_ratio > 0.99
        @test_throws ArgumentError map_gabls3_to_jacobian(1, [0.1], [0.1, 0.2], [0.01], config)
        @test_throws ArgumentError BifurcationGatingParams(epsilon_on=-0.05, epsilon_off=-0.10)
    end
end