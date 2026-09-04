using Test  
using LinearAlgebra  
using SBLGating  
  
@testset "SBLGating Synthetic Benchmarks" begin  
    params = BifurcationGatingParams(  
        epsilon_on = 1e-4,  
        epsilon_off = -1e-4,  
        zeta_z_tol = 1e-2,  
        km_floor = 0.1,  
        kh_floor = 0.01  
    )  
  
    @testset "Scenario 1: Coordinate Fold Passage (LLJ Nose)" begin  
        state = GatingState(Float64)  
        # Stable fast mode (λ_f > **ϵ**_on) with coordinate fold (|ζ_z| → 0)  
        J_fast = [5e-4 0.1; 0.0 -0.05]  
        lambda_f = extract_fast_eigenvalue(J_fast, state)  
        zeta_z = 1e-4  # Coordinate fold condition satisfied  
  
        Km, Kh, gated = update_gating_state!(state, lambda_f, zeta_z, 0.005, 0.001, params)  
  
        @test gated == true  
        @test Km == params.km_floor  
        @test Kh == params.kh_floor  
    end  
  
    @testset "Scenario 2: Genuine Dynamical Bifurcation" begin  
        state = GatingState(Float64)  
        state.is_active = true  # Gate initially active  
  
        # Collapse of normal hyperbolicity (λ_f < **ϵ**_off)  
        J_fast = [-5e-3 0.0; 0.0 -0.1]  
        lambda_f = extract_fast_eigenvalue(J_fast, state)  
        zeta_z = 1e-4  
  
        Km, Kh, gated = update_gating_state!(state, lambda_f, zeta_z, 0.005, 0.001, params)  
  
        @test gated == false  
        @test Km == 0.005  # Standard diagnostic diffusivity preserved  
        @test Kh == 0.001  
    end  
  
    @testset "Scenario 3: Hysteresis Chatter Suppression" begin  
        state = GatingState(Float64)  
  
        # 1. Trigger Gate ON  
        update_gating_state!(state, 2e-4, 0.0, 0.01, 0.001, params)  
        @test state.is_active == true  
  
        # 2. Fluctuate in deadband (**ϵ**_off < λ_f < **ϵ**_on) -> Gate remains ON  
        _, _, gated = update_gating_state!(state, 0.0, 0.0, 0.01, 0.001, params)  
        @test gated == true  
  
        # 3. Cross negative threshold (λ_f < **ϵ**_off) -> Gate turns OFF  
        _, _, gated = update_gating_state!(state, -2e-4, 0.0, 0.01, 0.001, params)  
        @test gated == false  
    end  
  
    @testset "Scenario 4: Eigenvector Continuation Across Crossings" begin  
        state = GatingState(Float64)  
          
        # Initial step: fast eigenvector along [1.0, 0.0]  
        J1 = [2e-4 0.0; 0.0 -1e-2]  
        λ1 = extract_fast_eigenvalue(J1, state)  
          
        # Step 2: Swap diagonal order in matrix; continuation must track fast mode [1.0, 0.0]  
        J2 = [-1e-2 0.0; 0.0 2e-4]  
        λ2 = extract_fast_eigenvalue(J2, state)  
  
        @test λ1 ≈ 2e-4  
        @test λ2 ≈ 2e-4  
    end  
end  
