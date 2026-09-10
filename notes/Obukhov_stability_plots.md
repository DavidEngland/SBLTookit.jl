The two Julia scripts—**`plot_obukhov_heatmaps.jl`** and **`plot_stability_heatmaps.jl`**—share a synchronized \(2 \times 2\) spatiotemporal diagnostic layout, but they differ fundamentally in their **underlying physical inputs, mathematical closure models, vertical domain extents, and numerical regularization techniques**.

Below is an explicit comparison of the two engines as implemented in the codebase and visualized in **`heatmaps.pdf`**:

---

### 1. Primary Physical Driving Input (Fluxes vs. State Gradients)

* **`plot_obukhov_heatmaps.jl` (Flux-Driven):**
  * **Input Mechanism:** Ingests surface turbulent flux measurements—specifically friction velocity (\(u_*\)) and surface heat flux (\(H_s\) or \(\overline{w'\theta'}\)).
  * **Stability Calculation:** Directly evaluates the Monin–Obukhov length scale \(L = -\frac{u_*^3 \theta_0}{\kappa g \overline{w'\theta'}}\) and computes the similarity coordinate \(\zeta(z,t) = z / L(z,t)\) directly from turbulent transport quantities.
* **`plot_stability_heatmaps.jl` (Gradient-Driven):**
  * **Input Mechanism:** Ingests profile state gradients—specifically vertical wind shear (\(S\)) and potential temperature stratification (\(\partial \theta_v / \partial z\)) extracted from tower soundings or Single-Column Model (SCM) fields.
  * **Stability Calculation:** Evaluates the local gradient Richardson number profile \(Ri_g(z,t) = \frac{g}{\theta_0} \frac{\partial \theta_v / \partial z}{S^2}\) and **inverts** \(Ri_g \to \zeta\) through a similarity closure model.

---

### 2. Mathematical Closures & Inversion Architecture

* **`plot_obukhov_heatmaps.jl` (Classical / Businger–Dyer Inversion):**
  * Inverts stability using rational Businger–Dyer mappings (\(\zeta = \frac{Ri_g}{1 - 5 Ri_g}\)).
  * **Limitation:** Requires hard numerical caps (e.g., \(Ri_g \le 0.19\)) near critical thresholds (\(Ri_c \approx 0.20\)), which flatlines \(\zeta\) and artificially truncates hyper-stable nocturnal structure aloft.
* **`plot_stability_heatmaps.jl` (Non-Singular GLGS / eMOST Inversion):**
  * Employs the SHEBA-calibrated **Grachev et al. (2007) / GLGS** flux-profile relationships (\(\phi_m(\zeta) = 1 + \frac{a_m \zeta}{(1 + b_m \zeta)^{2/3}}\), \(\phi_h(\zeta) = Pr_0 (1 + \frac{a_h \zeta}{1 + b_h \zeta})\)).
  * Uses a 1D Newton-Raphson solver (`zeta_from_rig_glgs`) to invert \(Ri_g \to \zeta\) continuously across hyper-stable regimes (\(0 \le \zeta < 100\)) without numerical clipping or pole singularities. Panel 2 is explicitly titled **`$\zeta(z,t)$ [GLGS]`**.

---

### 3. Vertical Domain Extent and Grid Resolution

* **`plot_obukhov_heatmaps.jl`:** Standardized on a vertical evaluation grid of **\(z \in [1.0, 100.0\text{ m}]\)** at 1-meter increments (`COMMON_Z_GRID = 1.0:1.0:100.0`). Designed primarily for surface-layer and tower-footprint scaling.
* **`plot_stability_heatmaps.jl`:** Expanded to **\(z \in [1.0, 200.0\text{ m}]\)** at 1-meter increments (`COMMON_Z_GRID = 1.0:1.0:200.0`). This extended domain captures elevated Low-Level Jet (LLJ) cores, jet-driven shear erosion, and nocturnal boundary layer capping inversions (e.g., across the 200m Cabauw mast in GABLS3).

---

### 4. Derivative Stencils & Track A Regularization

* **`plot_obukhov_heatmaps.jl`:** Originally interpolated raw profiles onto a uniform 1m grid *before* applying 3-point uniform finite-difference operators (`uniform_gradient_1d`, `uniform_hessian_1d`).
  * *Artifact Note:* Piecewise-linear interpolation prior to second differentiation forces \(\zeta_{zz} = 0\) between sensor levels, producing horizontal "banded" Hessian zeroing artifacts.
* **`plot_stability_heatmaps.jl` (Track A Compliant):** Implements **`assemble_derivatives_track_a`**.
  * Constructs Taylor-Vandermonde derivative matrices (\(\mathbf{D}_1, \mathbf{D}_2\)) directly on the **native non-uniform tower levels** in log-height space (\(\xi = \ln(z/z_0)\)) before projecting smooth derivative fields onto the evaluation grid.
  * Enforces the **Mandatory \(N_z \ge 3\) Stencil Gate**; sparse datasets failing this gate (e.g., 2-level SHEBA) trigger a `StencilCollapseException` and are safely routed to a Bulk Richardson (\(Ri_b\)) fallback engine (`assemble_bulk_fallback`).

---

### 5. Summary of the 4 Diagnostic Panels in `heatmaps.pdf`

Both scripts output a \(2 \times 2\) panel matrix per campaign:

| Panel | Metric | Mathematical Formulation | Visual Target & Diagnostic Role |
| :--- | :--- | :--- | :--- |
| **Panel 1** | **Reciprocal Obukhov Length** | \(1/L(z,t) \quad [\text{m}^{-1}]\) | Symmetric `:coolwarm` / `:bwr` colormap (\(\pm 0.05\text{ m}^{-1}\)). Distinguishes convective instability (\(1/L < 0\)) from stable stratification (\(1/L > 0\)). |
| **Panel 2** | **Stability Parameter** | \(\zeta(z,t) = z / L(z,t)\) | Symmetric-log transform (\(\text{sgn}(\zeta)\log_{10}(1 + \|\zeta\|)\)) on a `:puor` colormap (\(\pm 1.5\)). Resolves near-neutral surface layers alongside hyper-stable regimes aloft. |
| **Panel 3** | **Coordinate Jacobian** | \(\zeta_z = \partial \zeta / \partial z \quad [\text{m}^{-1}]\) | Centered stencil first derivative. Superimposes a **solid white contour line at \(\zeta_z = 0\)** to explicitly map spatial coordinate fold loci (e.g., LLJ nose axes). |
| **Panel 4** | **Coordinate Curvature** | \(\zeta_{zz} = \partial^2 \zeta / \partial z^2 \quad [\text{m}^{-2}]\) | Centered stencil second derivative (\(\pm 0.02\text{ m}^{-2}\)). Maps inversion capping heights, shear erosion boundaries, and near-surface inflection masking zones. |

---

### Multi-Campaign Behavior in `heatmaps.pdf`

In `heatmaps.pdf`, evaluating these two scripts across campaign archetypes illustrates key resolution thresholds:

* **CASES-99 (\(N_z = 7, z \le 55\text{ m}\)):** Cleanly displays the descending nocturnal Low-Level Jet nose. Panel 3's white contour (\(\zeta_z = 0\)) tracks the coordinate fold locus tracking the jet core.
* **GABLS3 (\(N_z = 38, z \le 200\text{ m}\)):** Captures the full 24-hour diurnal cycle over the Cabauw mast, showing daytime convection (\(\zeta < 0\)) transitioning into elevated nocturnal shear erosion.
* **SHEBA (\(N_z = 2, z \in \{2.5, 10\}\text{ m}\)):** Triggers the \(N_z \ge 3\) stencil gate (`Bypassed Guard`), routing to the 1D bulk fallback layer to avoid finite-difference under-determination.
