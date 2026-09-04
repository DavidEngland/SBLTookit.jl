Use this sequential prompt sequence in your VS Code AI chat box to move systematically from structural setup to full A/B benchmark execution.  
**Phase 1: Planning & Environment Architecture**  
```
I am implementing a geometry-aware gating system for SBLToolkit.jl to address false runaway decoupling in stable boundary layer modeling. 

Please inspect my current workspace repository structure and:
1. Identify the exact directory location where new modules should live (e.g., `src/SBLGating.jl` and `src/GABLS3Adapters.jl`).
2. Check `Project.toml` to ensure `LinearAlgebra` and `Test` are present in dependencies.
3. Draft a step-by-step file creation and integration plan that avoids breaking existing GABLS3 single-column model (SCM) runs.

```
**Phase 2: Core Module Implementation**  
```
Please create two new Julia files in `src/`:

1. `src/SBLGating.jl`:
   - Implement `BifurcationGatingParams` struct (with hysteresis thresholds `epsilon_on`, `epsilon_off`, and `zeta_z_tol`).
   - Implement `GatingState` mutable struct tracking `is_active` and `prev_v_f`.
   - Implement `extract_fast_eigenvalue` using eigenvector continuation (dot product overlap) to avoid mode crossing.
   - Implement `update_gating_state!` decoupling normal hyperbolicity checks from geometric regularity (`|zeta_z| <= zeta_z_tol`), applying mixing floors `km_floor` and `kh_floor` only when gated.

2. `src/GABLS3Adapters.jl`:
   - Implement `GSPTModelConfig` for timescale separation (`epsilon`), regularization scale (`epsilon_c`), and length scale (`l0`).
   - Implement `map_gabls3_to_jacobian` mapping SCM profiles (E, S, N^2) to a 2x2 scaled Jacobian J. Use C1-continuous Hill-type buoyancy flux derivatives.
   - Implement `compute_fold_ratio(zeta_z, zeta_zz, C_constitutive)` to compute mapping curvature dominance C_mapping > 0.99.

Ensure full type stability (`T<:AbstractFloat`) across all functions.

```
**Phase 3: Synthetic Unit Testing**  
```
Create a test file at `test/test_sbl_gating.jl` and include it in `test/runtests.jl`.

Implement `@testset` cases using parameters calibrated to CASES-99 Night 991018:
1. Fold Illusion (LLJ Nose): Set J_fast to yield a strongly attracting fast mode (lambda_f = -0.412 s⁻¹) and a coordinate profile knee (zeta_z ≈ 0). Assert `gated == true` and mixing floors are applied.
2. Genuine Bifurcation: Set lambda_f = -0.02 s⁻¹ (> epsilon_off). Assert `gated == false` and standard diagnostic diffusivity is preserved.
3. Deadband Hysteresis: Test transition across epsilon_on (-0.10) and epsilon_off (-0.05) to verify chatter suppression.
4. Mode Continuation: Test eigenvector tracking across diagonal swaps in J_fast.

Run the test suite and verify all assertions pass cleanly.

```
**Phase 4: SCM Driver Integration & A/B Validation**  
```
Now hook `SBLGating` and `GABLS3Adapters` into the main SCM driver in `SBLToolkit.jl`:

1. Implement `run_gabls3_step!` to iterate over height levels, construct J_fast, extract lambda_f, compute the Fold Ratio audit, and update gating states before updating vertical diffusion matrices.
2. Create an A/B benchmark script `scripts/run_gabls3_ab_experiment.jl`:
   - Experiment A: Baseline GABLS3 run (`enable_gating = false`).
   - Experiment B: Geometry-Aware Gated run (`enable_gating = true`).
3. Output comparative diagnostics: Plot/print the 24-hour surface temperature (T_2m) time series, boundary layer height, and total gated timesteps. Compute the exact reduction in the 3.5 K surface cold bias.

```
Should we start with the Phase 1 workspace audit prompt, or are your file paths already structured and ready for Phase 2?  
