using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
using SBLToolkit, SBLToolkit.GSPTPhase2

using Printf, Statistics, LinearAlgebra

function run_formal_phase2_benchmark()
        println("="^100)
        println(" SBLToolkit.jl: Executing Primitive-Space MDP GSPT Phase 2 Diagnostic Pipeline")
        println("="^100)

        # Cabauw 9-level mast configuration (m)
        z_grid = [2.0, 10.0, 20.0, 40.0, 80.0, 120.0, 140.0, 180.0, 200.0]
        n = length(z_grid)

        # Synthetic profile with Low-Level Jet (LLJ) nose at 120m (S^2 -> 0)
        u_base = [2.0, 5.0, 8.0, 12.0, 15.0, 15.01, 14.2, 11.0, 9.0] # LLJ nose around 120m
        v_base = [0.2, 0.5, 0.8, 1.2, 1.5, 1.50, 1.4, 1.1, 0.9]
        th_base = 285.0 .+ 0.04 .* z_grid
        wth_base = -0.015 .* exp.(-z_grid ./ 50.0)
        uw_base = -0.05 .* exp.(-z_grid ./ 60.0)
        vw_base = -0.01 .* exp.(-z_grid ./ 60.0)
        tke_base = 0.2 .* exp.(-z_grid ./ 40.0) .+ 0.01
        Km_base = 0.5 .* (z_grid ./ 10.0) .* exp.(-(z_grid ./ 30.0) .^ 2) .+ 0.001

        # Simulated observational noise (5 cm/s on wind, 0.05 K on temp)
        obs_data = ProfileData(z_grid, u_base, v_base, th_base, wth_base,
                uw_base, vw_base, tke_base, Km_base,
                0.05, 0.05, 0.05, 0.01)
        obs_res = GSPTPhase2.compute_gspt(obs_data; is_observation=true, S2_min=1e-3)

        # SCM profile (Over-smoothed wind gradient)
        scm_u = [2.0, 4.0, 6.0, 9.0, 12.0, 13.5, 14.0, 13.0, 11.0]
        scm_data = ProfileData(z_grid, scm_u, v_base, th_base, wth_base,
                uw_base, vw_base, tke_base, Km_base, 0.0, 0.0, 0.0, 0.0)
        scm_res = GSPTPhase2.compute_gspt(scm_data; is_observation=false, S2_min=1e-3)

        # LES profile
        les_data = ProfileData(z_grid, u_base, v_base, th_base, wth_base,
                uw_base, vw_base, tke_base, Km_base, 0.0, 0.0, 0.0, 0.0)
        les_res = GSPTPhase2.compute_gspt(les_data; is_observation=false, S2_min=1e-3)

        # Tangential Cone Condition Check across SBL perturbation
        x1 = [u_base; v_base; th_base]
        x2 = x1 .+ 0.02 .* randn(length(x1))
        D1_dim, _ = GSPTPhase2.build_operators(z_grid)
        gamma_val = GSPTPhase2.check_tangential_cone(x1, x2, 9.81, 285.0, D1_dim)

        # Track comparison
        metrics = GSPTPhase2.compare_tracks(obs_res, scm_res, les_res)

        # Print Detailed Diagnostics
        @printf("%-6s | %-12s | %-12s | %-12s | %-12s | %-12s\n",
                "z (m)", "R_coord(Obs)", "R_coord(SCM)", "B_R(SCM)", "Δ_close(Obs)", "τ_reg_sens")
        println("-"^95)
        for i in 1:n
                flag = obs_res.obs_diag.ill_conditioned_mask[i] ? "*" : " "
                @printf("%6.1f%s| %12.4f | %12.4f | %12.4f | %12.6f | %12.6f\n",
                        z_grid[i], flag, obs_res.const_geom.R_coord[i], scm_res.const_geom.R_coord[i],
                        metrics.bias_R_scm[i], obs_res.obs_diag.closure_residual[i], obs_res.obs_diag.tau_reg_sens[i])
        end
        println("="^95)
        println("(* indicates ill-conditioned level masked via S^2 <= S^2_min threshold)")
        @printf("Operator Grid Condition Number κ(R_tilde): %.2e\n", obs_res.obs_diag.grid_cond_number)
        @printf("Empirical Tangential Cone Residual Γ(x1, x2): %.4f\n", gamma_val)
        @printf("Inflection Masking Fraction (|R+1| < 0.15): Obs = %.2f%%, SCM = %.2f%%, LES = %.2f%%\n",
                metrics.masking_fraction_obs * 100, metrics.masking_fraction_scm * 100, metrics.masking_fraction_les * 100)
        @printf("SCM Coordinate Bias Metrics: RMSE = %.4f, MAE = %.4f, Median |ΔR| = %.4f\n",
                metrics.rmse_scm, metrics.mae_scm, metrics.median_abs_diff_scm)
        @printf("LES Coordinate Bias Metrics: RMSE = %.4f, MAE = %.4f, Median |ΔR| = %.4f\n",
                metrics.rmse_les, metrics.mae_les, metrics.median_abs_diff_les)
        println("="^95)
end

run_formal_phase2_benchmark()