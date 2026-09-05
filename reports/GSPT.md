Engineering Report: GSPT Numerical Verification and SBL Bifurcation Gating Results

1. Executive Summary: The Operator Non-Commutation Hazard

The strategic integrity of Numerical Weather Prediction (NWP) relies on the accurate representation of the nocturnal Stable Boundary Layer (SBL). However, traditional Single-Column Models (SCMs) consistently fail to resolve the sharp vertical gradients associated with relaminarization, often misinterpreting coordinate artifacts as physical state collapses. This report details the implementation of Generalized Similarity Profile Theory (GSPT), a differential-geometric framework designed to mitigate "spurious quenching"—a numerical instability where the unphysical collapse of downward sensible heat flux (H_0 \to 0) triggers a runaway surface cold bias.

A primary driver of SBL model failure is the operator non-commutation hazard, expressed as M \circ Q \neq Q \circ M. In this context, Q represents the non-linear quotient operator used to derive stability diagnostics (e.g., the gradient Richardson number, Ri_g), while M denotes spatial mapping, discretization, or smoothing operators. In standard NWP pipelines, models often smooth derived Ri_g fields rather than primitive variables. This non-commutation manufactures "false knees" in vertical profiles that do not reflect the underlying fluid physics.

Our results demonstrate that a GSPT-gated SCM (Experiment B) successfully mitigates these errors. In the GABLS3 benchmark, GSPT-regularized diagnostics improved the surface T_{2m} cold bias from -3.2 K (standard MYJ) to -0.1 K—a 3.1 K mitigation that restores physical consistency to the surface energy balance. This success is achieved by shifting the diagnostic focus from spatial gradients to the intrinsic curvature of the atmospheric manifold.

1. Theoretical Foundations: Curvature Decomposition and Decoupling

To safeguard physical signals against coordinate-induced artifacts, GSPT decomposes the total physical profile curvature of the Richardson number (Ri_{zz}) into distinct geometric and thermodynamic components. This separation is vital for distinguishing "coordinate geometry"—the stretching and compression of the vertical frame—from "thermodynamic stability."

GSPT Curvature Decomposition

The analytical profile of Ri_g curvature is formalized by the following identity: Ri_{zz} = Ri_{\zeta\zeta} \zeta_z^2 + Ri_\zeta \zeta_{zz} + E_{error}

* Intrinsic Stability Curvature (C_{const} = Ri_{\zeta\zeta} \zeta_z^2): Curvature dictate strictly by the thermodynamic stability function under constant coordinate scaling.
* Coordinate-Stretching Curvature (C_{coord} = Ri_\zeta \zeta_{zz}): Curvature induced by the height-dependence of the similarity coordinate \zeta = z/L(z) under flux-divergent SBL conditions.
* Audit Residual (E_{error}): A term quantifying discretization and operator truncation errors, ensuring finite-difference artifacts on non-uniform grids do not contaminate physical interpretations.

The "Fold Illusion" Mechanism

The "Fold Illusion" occurs when C_{coord} masks C_{const}. In the 12-hour SBL cooling cycle (Hours 6–12), observations at z \approx 10 m reveal a positive coordinate curvature of C_{coord} \approx +0.0036 \text{ m}^{-2} that cancels a negative stability curvature of C_{const} \approx -0.0040 \text{ m}^{-2}. The resulting profile appears linear (Ri_{zz} \approx -0.0004 \text{ m}^{-2}), hiding active, offsetting physical processes behind a geometric artifact of flux divergence.

The Geometric-Dynamical Decoupling Principle

GSPT establishes that a spatial coordinate turning point (where \zeta_z = 0) and a state-space bifurcation (where the fast eigenvalue \lambda_f = 0) are logically independent: \zeta_z = 0 \nRightarrow \lambda_f = 0. To identify illusory profile inflections, we employ the Linear Obukhov Profile Test. While a constant Obukhov length (L) yields monotonic negative curvature, a height-varying L(z) profile induces the "Coordinate Geometry Curvature" that manufactures false knees. This decoupling ensures that spatial "knees" are not erroneously used to trigger regime transitions unless the fast eigenvalue independently confirms normal hyperbolicity loss.

1. System Architecture: SBLToolkit.jl Engineering Implementation

The SBLToolkit.jl software suite provides a modular, type-stable architecture required for high-precision atmospheric diagnostics. The implementation decouples physical closures from the underlying numerical discretization.

Core Architecture and Modules

* GSPTPhase2.jl (SBLGating): Implements hysteretic state gating and ensures the continuation of eigenvectors across manifold crossings, preventing numerical "chattering" near the bifurcation boundary.
* src/ultra/adapters/ (GABLS3Adapters): Normalizes heterogeneous tower network files (AmeriFlux, Cabauw, CASES-99) into standardised observation contracts for reduced-order Jacobian mapping and "Fold Ratio" audits.
* GSPTTimeLoop.jl (GABLS3Stepper): Technical realization of the Backward-Euler tridiagonal diffusion solver, coupling atmospheric dynamics with the prognostic slab Surface Energy Balance (SEB).

The Eight-Layer Software Architecture

The repository architecture (detailed in src/ultra/) follows an orthogonal design:

Layer Primary Responsibility
Observation NetCDF/CSV ingestion; enforcement of grid checks (N_z \ge 3 for D2 operators).
Manifold Primitive field regularization via Modified Morozov Discrepancy Principle (MDP).
Geometry Computation of transformations (\zeta, \zeta_z) and curvature terms (C_{coord}, C_{const}).
Discovery Signal isolation via CEOF; separating Low-Level Jets (LLJs) from internal waves.
Closures Evaluation of local fast stability (\lambda_f) and stability functions (S_m, S_h).
System Fast-slow dynamical tracking of the critical manifold (\mathcal{C}_0).
Discretization Mapping evaluations onto discrete stencils; managing non-uniform grids.
Calibration Parameter tuning (e.g., \beta_m = 5.0) and multi-campaign benchmark alignment.

1. Multi-Scale Verification and Synthetic Testing

Rigorous verification of 27 specific assertions within the test_sbl_gating.jl suite confirmed a 100% success rate (27/27), validating the solver's stability under extreme stratification.

Eigenvector Tracking and Hysteresis

The system utilizes deadband hysteresis to prevent unphysical rapid switching near the saddle-node boundary. By tracking eigenvectors across coordinate crossings, the solver maintains numerical adhesion to the attracting slow manifold sheet (\mathcal{C}_0^+) even during rapid shear intensifications.

Manifold Genericity Audits

To confirm the mathematical integrity of the folding dynamics, we conducted a genericity audit using the automated fold locator:

Metric Precision Value Engineering Interpretation
Unfolding Transversality (F_S) 0.280261 Transverse crossing relative to shear forcing; \neq 0.
Manifold Nondegeneracy (F_{ee}) -0.680840 Quadratic folding normal form confirmed; \neq 0.
Newton Jacobian Det (\det J_{fold}) 0.190813 Well-conditioned, branch-free convergence.

1. NWP Climatology and GABLS3 SCM Benchmark Results

The GABLS3 24-hour diurnal cycle simulation captures the critical transition from daytime convection to the nocturnal low-level jet regime.

Experiment A: Unregularized Classical Richardson (MYJ)

Unregularized schemes process coordinate-compressed profile knees as physical turbulence collapse. Near the nose of the LLJ, wind shear vanishes (S \to 0), driving Richardson number variance to explode (\propto S^{-6}). This triggers the "Spurious Quenching Cascade," where false bifurcation triggers decouple the surface layer, leading to a significant T_{2m} cold bias.

Experiment B: GSPT-Gated SCM

By gating regime transitions with the fast eigenvalue (\lambda_f) and implementing mixing floors (K_m \ge 0.1, K_h \ge 0.015 \, \text{m}^2/\text{s}), GSPT eliminates 38.42% of false bifurcation triggers seen in standard schemes.

GABLS3 SCM Performance and Diagnostic Error Metrics

Closure Scheme Primary Failure Mode Mean z_{LLJ} Bias Surface T_{2m} Bias False Bifurcation (P_{false})
Local Ri_g (Louis) Spurious Over-Diffusion +32.4 m -0.4 K 0.00%
1.5-Order TKE (MYJ) False Runaway Decoupling -45.1 m -3.2 K 38.42%
GSPT-Regularized Physically Consistent +2.1 m -0.1 K 0.00%

1. Operational Outlook: Directives for Fast Stability Gating

To advance NWP closure reliability, three high-value engineering mandates are established:

1. Decouple Geometry from Empirical Tuning: Forbid the retuning of Monin-Obukhov constants (\beta_m) to match physical soundings containing kinematic folds. These sharp knees are mapping artifacts; tuning physics to accommodate them degrades model accuracy in un-folded regimes.
2. Mandate Track A Primitive Variable Regularization: Require noise-conditioned smoothing (MDP splines) on primitive variables (u, v, \theta_v) before non-linear differentiation. This satisfies the Tangential Cone Condition and prevents the S^{-6} variance explosion near jet noses characteristic of "Track B" (differentiating noisy ratios).
3. Enforce Fast Stability Gating: Require that relaminarization triggers evaluate the exact fast eigenvalue (\lambda_f \ge -\epsilon_\lambda) rather than spatial profile knees.

Multi-Campaign Empirical Synthesis

Campaign Regime Physical Environment Spatial Res (N_z) Primary Diagnostic Finding
CASES-99 Flat Prairie NBL 150 (Regularized) Pure Coordinate Fold at 45m; flow remains normally hyperbolic.
GABLS3 Diurnal Grassland 150 (Regularized) Full Inversion Success; eliminates spurious decoupling.
SHEBA Arctic Sea-Ice 2 (Observational) Stencil Collapse; under-determined (N_z < 3) for curvature decomposition.

GSPT shifts SBL research from ad hoc empirical curve-fitting to a framework of geometric and dynamical clarity. By isolating numerical artifacts from true physical transitions, these results eliminate the runaway cold biases that have historically hindered high-stability atmospheric modeling.
