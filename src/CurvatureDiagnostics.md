We will proceed with **Option A: Stage 4 Ingestion and GCV Spline Loop**, as clean primitive gradient extraction is the necessary prerequisite before tracking phase-space trajectories in Option B.

Below is an operational Julia module for `SBLToolkit.jl` that enforces operator non-commutation by applying Generalized Cross-Validation (GCV) spline smoothing directly to primitive fields ($u, v, \theta_v$) prior to non-linear combination and differentiation.

**Key Features of the Code**

* **Primitive Smoothing Operator:** Fits `SmoothingSpline` stencils directly to raw $u, v, \theta_v$ arrays before calculating shears or buoyancy gradients, mathematically preventing $\mathcal{O}(U_z^{-6})$ noise blowup near low-level jet noses.
* **Inverse Obukhov Coordinate:** Builds the finite inverse-length profile $\chi(z) = 1/L(z)$, smooths $\chi$ directly, and obtains $\zeta_z = \chi + z\chi'$ and $\zeta_{zz} = 2\chi' + z\chi''$ without evaluating a divergent $L$ near neutral conditions.
* **Direct Flux Override:** When `tau_raw` and `heat_flux_raw` are supplied, uses $\chi(z) = -\kappa(g/\theta_0)\overline{w'\theta_v'}(z)/\tau(z)^{3/2}$ instead of the gradient Richardson proxy.
* **CASES-99 CSV Adapter:** `process_curvature_csv("workspace/out/gspt_cases99_coordinates.csv")` groups the flat trajectory by timestamp, treats its `Ri_g` column as the documented $\chi$ proxy, preserves its supplied `C_M`, and returns time-by-height diagnostic matrices.
* **Domain-Adaptive Regularization:** Evaluates $K_0$ as $\text{median}(\vert{}Ri_{g,zz}\vert{}) + 10^{-6}$ per profile, preserving dimensionally consistent scaling $[K_0] = \text{m}^{-2}$ across varying stability regimes.
* **Automated State Flagging:** Maps output points into the five-state taxonomy (`:PureCoordinateFold`, `:DynamicSingularityCandidate`, `:VerifiedSaddleNode`, `:HybridFold`, `:ClosureBreakdown`). Optional `state_singularity`, `verified_saddle_node`, and `closure_residual` inputs provide independent state-space and closure evidence.

Are your DNS data arrays structured as binary flat files or standard NetCDF/HDF5, and should we integrate an ingestion interface for those file formats next?
