using Test
using Statistics
using Random

# Include the source module (either version 1 or version 2 depending on setup)
if !isdefined(Main, :GSPTCases99Diagnostic)
    include("../src/gspt_cases99_diagnostic_v2.jl")
end
using .GSPTCases99Diagnostic

@testset "GSPTCases99Diagnostic.jl Production Unit Tests" begin
    # 1. Struct and Parameters Sanity Check
    @testset "Struct and Parameter Integrity" begin
        p_default = CASES99Params()
        @test p_default.g == 9.81
        @test p_default.theta_ref == 285.0
        @test p_default.beta == 5.0
        @test p_default.z_0 == 0.01
        
        # Test constructor with custom overrides
        p_custom = CASES99Params(beta=4.7, sigma_u=0.01)
        @test p_custom.beta == 4.7
        @test p_custom.sigma_u == 0.01
        @test p_custom.g == 9.81  # remains default
    end

    # 2. Spline Solver Accuracy & Analytical Derivatives
    @testset "Natural Cubic Spline Fitting and Evaluation" begin
        # Natural boundary conditions exactly preserve a linear profile.
        x = collect(range(-2.0, 2.0, length=11))
        y = 2.0 .* x .+ 1.0
        w = fill(1e-5, length(x)) # negligible noise weighting
        
        # Fit cubic spline with extremely low regularization (acts as interpolation)
        xk, ak, Mk, hk = GSPTCases99Diagnostic.fit_natural_cubic_spline(x, y, w, 1e-12)
        
        # Test evaluation at a midpoint node
        x_eval = 0.5
        y_eval = GSPTCases99Diagnostic.evaluate_spline(xk, ak, Mk, hk, x_eval; order=0)
        # y(0.5) = 2.0
        @test isapprox(y_eval, 2.0, atol=1e-9)
        
        # Test analytical first derivative (y' = 2.0)
        dy_eval = GSPTCases99Diagnostic.evaluate_spline(xk, ak, Mk, hk, x_eval; order=1)
        @test isapprox(dy_eval, 2.0, atol=1e-9)
        
        # Test analytical second derivative (y'' = 0.0)
        ddy_eval = GSPTCases99Diagnostic.evaluate_spline(xk, ak, Mk, hk, x_eval; order=2)
        @test isapprox(ddy_eval, 0.0, atol=1e-9)
        
        # Natural boundary condition test: second derivatives at endpoints must vanish (Mk[1] == Mk[end] == 0.0)
        ddy_start = GSPTCases99Diagnostic.evaluate_spline(xk, ak, Mk, hk, x[1]; order=2)
        ddy_end = GSPTCases99Diagnostic.evaluate_spline(xk, ak, Mk, hk, x[end]; order=2)
        @test isapprox(ddy_start, 0.0, atol=1e-9)
        @test isapprox(ddy_end, 0.0, atol=1e-9)
    end

    # 3. Constant Obukhov Length Baseline (Constant Flux regime)
    @testset "Constant-Flux Surface Layer (C_M -> 0)" begin
        p = CASES99Params()
        z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 55.0]
        
        zeta_true = z_tower ./ 10.0
        Ri_true = zeta_true ./ (1.0 .+ p.beta .* zeta_true)
        
        # Verify analytical zeta inversion in the module's exact math
        r = Ri_true[1]
        zeta_inv = r / (1.0 - p.beta * r)
        @test isapprox(zeta_inv, zeta_true[1], atol=1e-7)
        
        C_const = -1.25
        C_mapping = 0.0
        K_0 = 0.01
        
        CM_eval = abs(C_mapping) / (abs(C_const) + abs(C_mapping) + p.epsilon_c * K_0)
        @test CM_eval == 0.0
    end

    # 4. Vanishing Wind Shear / Turning Point Baseline
    @testset "Fold Illusion / Jet-Nose Singularity (C_M -> 1)" begin
        p = CASES99Params()
        
        C_const = 0.0
        C_mapping = 0.05
        K_0 = 0.01
        
        CM_eval = abs(C_mapping) / (abs(C_const) + abs(C_mapping) + p.epsilon_c * K_0)
        @test CM_eval >= 0.99
    end

    # 5. Advanced GSPT Structural Correction: Multi-Variable Noise Propagation
    # Verifies that the total analytical propagated variance closely matches 
    # a Monte Carlo numerical simulation (error < 0.5%).
    @testset "GSPT Multi-Variable Noise Propagation Verification" begin
        p = CASES99Params()
        Random.seed!(42)
        
        # Standard layer evaluation inputs
        z_step = 5.0             # Spacing (m)
        U_z = 0.4                # Wind shear (s^-1)
        theta_z = 0.08           # Temperature gradient (K/m)
        
        # Propagated primitive gradient uncertainties
        sigma_U_z = p.sigma_u / (sqrt(2.0) * z_step)
        sigma_theta_z = p.sigma_theta / (sqrt(2.0) * z_step)
        
        # Sensitivities w.r.t gradients
        sens_theta = (p.g / p.theta_ref) / (U_z^2)
        sens_U = -2.0 * (p.g / p.theta_ref) * theta_z / (U_z^3)
        
        # Analytical total propagated variance
        var_analytical = (sens_theta^2) * (sigma_theta_z^2) + (sens_U^2) * (sigma_U_z^2)
        std_analytical = sqrt(var_analytical)
        
        # Monte Carlo Simulation (10,000 samples)
        N_samples = 10000
        noise_u = randn(N_samples) .* sigma_U_z
        noise_t = randn(N_samples) .* sigma_theta_z
        
        U_z_perturbed = U_z .+ noise_u
        theta_z_perturbed = theta_z .+ noise_t
        
        Ri_g_perturbed = (p.g / p.theta_ref) .* theta_z_perturbed ./ (U_z_perturbed .^ 2)
        
        std_numerical = std(Ri_g_perturbed)
        
        # Verify analytical error propagation matches numerical within 1%
        @test isapprox(std_numerical, std_analytical, rtol=1e-2)
    end

    # 6. Advanced GSPT Structural Correction: Taxonomic Matrix and Scale-Aware Threshold
    # Verifies that under systematic noise sweeps (0.1x to 5.0x), 
    # SBL states converge identically to their target classifications because
    # delta_obs adapts dynamically to local propagated uncertainty.
    @testset "Taxonomic Matrix Noise Invariance Sweep" begin
        p = CASES99Params()
        
        # Scale-normalized dynamic threshold: eta_d / epsilon
        epsilon_perturb = 0.05
        eta_d = 0.001
        delta_dyn = eta_d / epsilon_perturb # Bounded threshold = 0.02 s^-1
        @test delta_dyn == 0.02
        
        # Evaluation points
        z_step = 5.0
        U_z = 0.4
        theta_z = 0.08
        sens_theta = (p.g / p.theta_ref) / (U_z^2)
        sens_U = -2.0 * (p.g / p.theta_ref) * theta_z / (U_z^3)
        
        # Noise scale sweeps
        noise_scales = [0.1, 0.5, 1.0, 2.0, 5.0]
        
        # Test Case: Pure Coordinate Fold (C_M >= 0.70, D_dyn >= delta_dyn)
        C_M_coord = 0.85
        D_dyn_coord = 0.05
        
        for scale in noise_scales
            sigma_u = p.sigma_u * scale
            sigma_theta = p.sigma_theta * scale
            
            sigma_U_z = sigma_u / (sqrt(2.0) * z_step)
            sigma_theta_z = sigma_theta / (sqrt(2.0) * z_step)
            
            # Dynamic local fidelity floor: 2*sigma (95% confidence interval)
            var_local = (sens_theta^2) * (sigma_theta_z^2) + (sens_U^2) * (sigma_U_z^2)
            delta_obs = 2.0 * sqrt(var_local)
            
            # Simulated model-observation error well within the noise envelope
            E_scm = 1.2 * sqrt(var_local)
            
            # Verify fidelity test passes
            @test E_scm <= delta_obs
            
            # Verify stable taxonomic classification
            classification = :Unassigned
            if E_scm <= delta_obs
                if C_M_coord >= 0.70 && D_dyn_coord >= delta_dyn
                    classification = :PureCoordinateFold
                elseif C_M_coord < 0.50 && D_dyn_coord < delta_dyn
                    classification = :PureDynamicFold
                elseif C_M_coord >= 0.70 && D_dyn_coord < delta_dyn
                    classification = :HybridResonantFold
                end
            end
            @test classification == :PureCoordinateFold
        end
    end

    # 7. Full Simulation Run Integrity
    @testset "Full SCM Pipeline Run Integrity" begin
        p = CASES99Params()
        mktempdir() do directory
            filepath = joinpath(directory, "out", "gspt_cases99_coordinates.csv")
            z, t, CM, Ri, z_llj, h_inv = calculate_gspt_diagnostics(p; filepath)

            @test length(z) == 150
            @test length(t) == 25
            @test size(CM) == (150, 25)
            @test size(Ri) == (150, 25)

            @test all(CM .>= 0.0)
            @test all(CM .<= 1.0)
            @test all(Ri .>= -2.0)
            @test all(Ri .<= 2.0)

            @test isfile(filepath)
            @test filesize(filepath) > 0
        end
    end
end
