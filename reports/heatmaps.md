### 1. Primary Physical Driving Input (Fluxes vs. State Gradients)

* **`plot_obukhov_heatmaps.jl` (Flux-Driven):**
* **Input Mechanism:** Ingests surface turbulent flux measurements—specifically friction velocity ($u_*$) and surface heat flux ($H_s$ or $\overline{w'\theta'}$).
* **Stability Calculation:** Directly evaluates the Monin–Obukhov length scale $L = -\frac{u_*^3 \theta_0}{\kappa g \overline{w'\theta'}}$ and computes the similarity coordinate $\zeta(z,t) = z / L(z,t)$ directly from turbulent transport quantities.

* **`plot_stability_heatmaps.jl` (Gradient-Driven):**
* **Input Mechanism:** Ingests profile state gradients—specifically vertical wind shear ($S$) and potential temperature stratification ($\partial \theta_v / \partial z$) extracted from tower soundings or Single-Column Model (SCM) fields.
* **Stability Calculation:** Evaluates the local gradient Richardson number profile $Ri_g(z,t) = \frac{g}{\theta_0} \frac{\partial \theta_v / \partial z}{S^2}$ and **inverts** $Ri_g \to \zeta$ through a similarity closure model.

---

### 2. Mathematical Closures & Inversion Architecture

* **`plot_obukhov_heatmaps.jl` (Direct / Classical Businger–Dyer):**
* Evaluates $\zeta = z/L$ directly from surface fluxes without an inversion step.
* When profile inversions are evaluated, rational Businger–Dyer mappings ($\zeta = \frac{Ri_g}{1 - 5 Ri_g}$) require hard numerical caps (e.g., $Ri_g \le 0.19$) near critical thresholds ($Ri_c \approx 0.20$), flatlining $\zeta$ and artificially truncating hyper-stable nocturnal structure aloft.
* Panel 2 carries the descriptor **`$\zeta(z,t)$ [Direct]`** (or un-tagged) to signify pure flux-based scaling rather than a model inversion.

* **`plot_stability_heatmaps.jl` (Non-Singular GLGS / eMOST Inversion):**
* Employs SHEBA-calibrated **Grachev et al. (2007) / GLGS** flux-profile relationships:

$$\phi_m(\zeta) = 1 + \frac{a_m \zeta}{(1 + b_m \zeta)^{2/3}}, \quad \phi_h(\zeta) = Pr_0 \left(1 + \frac{a_h \zeta}{1 + b_h \zeta}\right)$$

* Uses a 1D Newton-Raphson solver (`zeta_from_rig_glgs`) to invert $Ri_g \to \zeta$ continuously across hyper-stable regimes ($0 \le \zeta < 100$) without numerical clipping or pole singularities. Panel 2 is explicitly titled **`$\zeta(z,t)$ [GLGS]`** to denote closure model inversion.
* **Visual Saturation Caveat:** While the Newton-Raphson solver converges up to $\zeta < 100$ without numerical clipping, the colormap scaling $\text{sgn}(\zeta)\log_{10}(1+\vert{}\zeta\vert{})$ capped at $\pm 1.5$ saturates visually at $\zeta \approx 30.62$ ($\log_{10}(1+30.62) = 1.5$). Calculated values in $30.6 \le \zeta \le 100$ render as saturated color extrema, creating visual truncation aloft despite continuous numerical output.

---

### 3. Vertical Domain Extent and Grid Resolution

* **`plot_obukhov_heatmaps.jl`:** Standardized on a vertical evaluation grid of **$z \in [1.0, 100.0\text{ m}]$** at 1-meter increments (`COMMON_Z_GRID = 1.0:1.0:100.0`). Designed primarily for surface-layer and tower-footprint scaling.
* **`plot_stability_heatmaps.jl`:** Expanded to **$z \in [1.0, 200.0\text{ m}]$** at 1-meter increments (`COMMON_Z_GRID = 1.0:1.0:200.0`). Captures elevated Low-Level Jet (LLJ) cores, jet-driven shear erosion, and nocturnal boundary layer capping inversions (e.g., across the 200m Cabauw mast in GABLS3).

---

### 4. Derivative Stencils & Track A Regularization

* **`plot_obukhov_heatmaps.jl`:** Interpolates raw profiles onto a uniform 1m grid *before* applying 3-point uniform finite-difference operators (`uniform_gradient_1d`, `uniform_hessian_1d`).
* *Artifact Note:* Piecewise-linear interpolation prior to second differentiation forces $\zeta_{zz} = 0$ between sensor levels, producing horizontal "banded" Hessian zeroing artifacts.

* **`plot_stability_heatmaps.jl` (Track A Compliant):** Implements **`assemble_derivatives_track_a`**.
* Constructs Taylor-Vandermonde derivative matrices ($\mathbf{D}_1, \mathbf{D}_2$) directly on the **native non-uniform tower levels** in log-height space ($\xi = \ln(z/z_0)$) before projecting smooth derivative fields onto the evaluation grid.
* Enforces the **Mandatory $N_z \ge 3$ Stencil Gate**; sparse datasets failing this gate (e.g., 2-level SHEBA) trigger a `StencilCollapseException` and are safely routed to a Bulk Richardson ($Ri_b$) fallback engine (`assemble_bulk_fallback`).

---

### 5. Summary of the 4 Diagnostic Panels

Both scripts output a $2 \times 2$ panel matrix per campaign:

| Panel | Metric | Mathematical Formulation | Visual Target & Diagnostic Role |
| --- | --- | --- | --- |
| **Panel 1** | **Reciprocal Obukhov Length** | $1/L(z,t) \quad [\text{m}^{-1}]$ | Symmetric `:coolwarm` / `:bwr` colormap ($\pm 0.05\text{ m}^{-1}$). Distinguishes convective instability ($1/L < 0$) from stable stratification ($1/L > 0$). |
| **Panel 2** | **Stability Parameter** | $\zeta(z,t) = z / L(z,t)$ | Symmetric-log transform ($\text{sgn}(\zeta)\log_{10}(1 + \Vert{}\zeta\Vert{})$) on a `:puor` colormap ($\pm 1.5$). Titled `[Direct]` for flux data and `[GLGS]` for inversions. Visually saturates above $\zeta \approx 30.6$. |
| **Panel 3** | **Coordinate Jacobian** | $\zeta_z = \partial \zeta / \partial z \quad [\text{m}^{-1}]$ | Centered stencil first derivative. Superimposes a **solid white contour line at $\zeta_z = 0$** to explicitly map spatial coordinate fold loci (e.g., LLJ nose axes). |
| **Panel 4** | **Coordinate Curvature** | $\zeta_{zz} = \partial^2 \zeta / \partial z^2 \quad [\text{m}^{-2}]$ | Centered stencil second derivative ($\pm 0.02\text{ m}^{-2}$). Maps inversion capping heights, shear erosion boundaries, and near-surface inflection masking zones. |

---

### Multi-Campaign Behavior

Evaluating these two scripts across campaign archetypes illustrates key resolution thresholds:

* **CASES-99 ($N_z = 7, z \le 55\text{ m}$):** Cleanly displays the descending nocturnal Low-Level Jet nose. Panel 3's white contour ($\zeta_z = 0$) tracks the coordinate fold locus tracking the jet core.
* **GABLS3 ($N_z = 38, z \le 200\text{ m}$):** Captures the full 24-hour diurnal cycle over the Cabauw mast, showing daytime convection ($\zeta < 0$) transitioning into elevated nocturnal shear erosion.
* **SHEBA ($N_z = 2, z \in \{2.5, 10\}\text{ m}$):** Triggers the $N_z \ge 3$ stencil gate (`Bypassed Guard`), routing to the 1D bulk fallback layer to avoid finite-difference under-determination.

---

## BLLAST (Flux-Driven Obukhov Diagnostic)

![BLLAST](./generated/sbltoolkit_heatmaps/bllast_obukhov_heatmaps.png)

## FLOSS II (Flux-Driven Obukhov Diagnostic)

![FLOSS II](./generated/sbltoolkit_heatmaps/floss_ii_obukhov_heatmaps.png)

## SHEBA (Flux-Driven Obukhov Diagnostic)

![SHEBA](./generated/sbltoolkit_heatmaps/sheba_obukhov_heatmaps.png)

## CASES-99 (Flux-Driven Obukhov Diagnostic)

![CASES-99](./generated/sbltoolkit_heatmaps/cases_99_obukhov_heatmaps.png)

## GABLS3 (Flux-Driven Obukhov Diagnostic)

![GABLS3](./generated/sbltoolkit_heatmaps/gabls3_obukhov_heatmaps.png)

## Code Reference: Obukhov Heatmaps

![CODE](../scripts/plot_obukhov_heatmaps.jl)
---

### Refined Architectural Summary

**1. Standardized Obukhov Heatmap Engine (`plot_obukhov_heatmaps.jl`)**
Interpolates heterogeneous tower observations onto a standard grid ($z \in [1, 100] \text{ m}$) at 1 m resolution (`COMMON_Z_GRID`).

* **$1/L(z,t)$ Reciprocal Length**: Bounded to $[-0.05, 0.05] \text{ m}^{-1}$ using a symmetric red-white-blue colormap (`:bwr`) to separate convective instability ($1/L < 0$) from stable boundary layers ($1/L > 0$).
* **$\zeta(z,t)$ Symmetric-Log Stability**: Employs $\text{sgn}(\zeta) \log_{10}(1 + \vert{}\zeta\vert{})$ on a Purple-White-Orange palette (`:puor`) across $[-1.5, 1.5]$ to resolve surface layer gradients alongside hyper-stable aloft conditions.
* **Coordinate Jacobian ($\zeta_z = \partial\zeta/\partial z$)**: Superimposes a solid white contour overlay at $\zeta_z = 0$ to identify coordinate fold locations.
* **Coordinate Curvature ($\zeta_{zz} = \partial^2\zeta/\partial z^2$)**: Maps vertical inversion capping layers and shear boundaries.

**2. GSPT Curvature Decomposition Engine (`gspt_curvature_contour_v3.jl`)**
Partitions total observed Richardson curvature $Ri_{zz}$ into constitutive thermodynamics and spatial coordinate geometry:

$$Ri_{zz} = \underbrace{R''(\zeta) \zeta_z^2}_{C_{\text{const}}} + \underbrace{R'(\zeta) \zeta_{zz}}_{C_{\text{coord}}} + \mathcal{E}_{\Delta z}$$

* **Non-Uniform Vandermonde Stencils ($\mathbf{D}_1, \mathbf{D}_2$)**: Solves local Taylor-Vandermonde linear systems ($\mathbf{A w} = \mathbf{b}$) to construct exact derivative operators for irregularly spaced tower levels ($z \in [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0] \text{ m}$).
* **Tikhonov Pre-Conditioning**: Filters high-frequency observational sensor noise ($\sigma_{\text{obs}} \approx 0.008$) prior to second differentiation via $(\mathbf{I} + \lambda_{\text{reg}} \mathbf{D}_2^T \mathbf{D}_2) Ri_{\text{smooth}} = Ri_{\text{obs}}$, ensuring the residual $\mathcal{E}_{\Delta z}$ reflects spatial discretization error rather than measurement noise.

**3. Primary Physical Diagnostic Takeaways**

* **Inflection Masking ($z \approx 10 \text{ m}$)**: Positive coordinate curvature ($C_{\text{coord}} \approx +0.0036 \text{ m}^{-2}$) cancels negative stability curvature ($C_{\text{const}} \approx -0.0040 \text{ m}^{-2}$), hiding active, opposing physical mechanisms behind an apparently linear profile ($Ri_{zz} \approx 0$).
* **Jet-Nose Fold Illusion ($z \approx 45 \text{ m}$)**: Near Low-Level Jet (LLJ) wind maxima where vertical shear vanishes ($S^2 \to 0$), visual profile knees are driven over 99.4% by coordinate-stretching geometry ($C_{\text{coord}}$) rather than true turbulence collapse.

---

## BLLAST (Gradient-Driven Stability Inversion)

![BLLAST](./generated/sbltoolkit_heatmaps/bllast_stability_heatmaps.png)

## FLOSS II (Gradient-Driven Stability Inversion)

![FLOSS II](./generated/sbltoolkit_heatmaps/floss_ii_stability_heatmaps.png)

## SHEBA (Gradient-Driven Stability Inversion)

![SHEBA](./generated/sbltoolkit_heatmaps/sheba_stability_heatmaps.png)

## CASES-99 (Gradient-Driven Stability Inversion)

![CASES-99](./generated/sbltoolkit_heatmaps/cases_99_stability_heatmaps.png)

## GABLS3 (Gradient-Driven Stability Inversion)

![GABLS3](./generated/sbltoolkit_heatmaps/gabls3_stability_heatmaps.png)

## Code Reference: Stability Heatmaps

![Code Reference](../scripts/plot_stability_heatmaps.jl)
