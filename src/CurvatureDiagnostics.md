We will proceed with **Option A: Stage 4 Ingestion and GCV Spline Loop**, as clean primitive gradient extraction is the necessary prerequisite before tracking phase-space trajectories in Option B.

Below is an operational Julia module for `SBLToolkit.jl` that enforces operator non-commutation by applying Generalized Cross-Validation (GCV) spline smoothing directly to primitive fields ($u, v, \theta_v$) prior to non-linear combination and differentiation.

**Key Features of the Code**

* **Primitive Smoothing Operator:** Fits `SmoothingSpline` stencils directly to raw $u, v, \theta_v$ arrays before calculating shears or buoyancy gradients, mathematically preventing $\mathcal{O}(U_z^{-6})$ noise blowup near low-level jet noses.
* **Analytical Spline Differentiation:** Evaluates $u_z, v_z, \theta_{v,z}, \zeta_z,$ and $\zeta_{zz}$ using spline basis derivatives (`predict(spl, z, order)`), ensuring $C^2$ continuity without introducing finite-difference truncation error.
* **Domain-Adaptive Regularization:** Evaluates $K_0$ as $\text{median}(\vert{}Ri_{g,zz}\vert{}) + 10^{-6}$ per profile, preserving dimensionally consistent scaling $[K_0] = \text{m}^{-2}$ across varying stability regimes.
* **Automated State Flagging:** Maps output points into the classification categories (`:PureCoordinateFold`, `:PureDynamicFold`, `:Ambiguous`, `:HybridResonantFold`).

Are your DNS data arrays structured as binary flat files or standard NetCDF/HDF5, and should we integrate an ingestion interface for those file formats next?