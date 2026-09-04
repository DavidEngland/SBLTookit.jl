# Nonlinear Closure–Discretization Instability and Spurious Quenching

Spurious turbulent quenching arises when spatially varying stability diagnostics interact with sharply nonlinear stability functions on a finite vertical grid. The underlying coordinate fold is a physically legitimate geometric feature; however, when its spatial curvature is under-resolved by the vertical grid, the discrete estimate of $Ri_g$ can cross a sharp stability threshold and trigger an unphysical numerical feedback loop.

---

### Mechanism Taxonomy

To prevent misattributing numerical artifacts to physical fold geometry, stability diagnostic behaviors are separated into three distinct pathways:

* **Physical Coordinate-Fold Curvature ($\zeta_z = 0, \, \zeta_{zz} \neq 0$):** A genuine geometric feature where $Ri_{g,zz}(z^*) = -R'(\zeta^*)\frac{z^* L''}{L^2}$. The physical field $Ri_g(z)$ remains continuous and smooth.
* **Discretization-Induced Threshold Crossing:** The physical fold curvature is legitimate, but the vertical grid resolution $\Delta z$ is insufficient to resolve $\frac{1}{2} Ri_{g,zz} (\Delta z)^2$. The discrete $Ri_g$ value crosses $Ri_{\text{crit}}$, driving an excessively large variation in eddy diffusivity $K_h$.
* **Gradient-Estimation Noise:** Differentiation of noisy flux estimates or wind/temperature shear fields ($Ri_g = N^2/S^2$) introduces numerical errors into $L(z)$ and its derivatives, artificially inflating local Richardson-number diagnostics prior to evaluation by the closure.

---

### Quadratic Sampling and Nonlinear Feedback Gain

Near a nondegenerate coordinate fold $z^*$ where $\zeta_z(z^*) = 0$ and $Ri_{g,z}(z^*) = 0$, a vertical grid centered on the fold samples the physical field as a locally quadratic profile:

$$Ri_g(z^* \pm \Delta z) = Ri_g(z^*) + \frac{1}{2} Ri_{g,zz}(z^*) (\Delta z)^2 + \mathcal{O}(\Delta z^3)$$

The fold generates a localized quadratic variation in $Ri_g$ rather than a mathematical discontinuity. Pathology arises when this variation interacts with the local linearization of the mixing scheme:

$$\delta K_h \simeq \mathcal{K}_h'(Ri_g) \, \delta Ri_g$$

In regions where $\vert{}\mathcal{K}_h'(Ri_g)\vert{} \gg 1$, a modest discrete variation $\delta Ri_g$ causes a disproportionate reduction in sensible heat flux $H = -\rho c_p K_h \theta_z$. This reduction isolates the lowest model layer, accelerating longwave radiative cooling and steepening the local potential temperature gradient $\theta_z$. Because $Ri_g \propto \theta_z$, the state perturbation feeds back directly into the subsequent stability diagnostic:

$$\delta Ri_g \longrightarrow \delta K_h \longrightarrow \delta H \longrightarrow \delta \theta_z \longrightarrow \delta Ri_g$$

The non-dimensional discrete feedback loop gain $\mathcal{G}$ is expressed as:

$$\mathcal{G} \sim \left(\frac{\partial Ri_g}{\partial \theta_z}\right) \left(\frac{\partial \theta_z}{\partial H}\right) \left(\frac{\partial H}{\partial K_h}\right) \left(\frac{\partial K_h}{\partial Ri_g}\right)$$

Spurious runaway quenching occurs when under-resolved grid sampling forces $\vert{}\mathcal{G}\vert{} \ge 1$, locking the lowest model layers into continuous, unphysical cooling.

---

### Operational Mitigation Strategies

Operational implementations must retain physical-space transport equations while regularizing the diagnostic evaluation of stability functions:

* **$C^2$-Consistent Profile Reconstruction:** Avoid naive finite-difference derivatives or low-pass spatial filters (e.g., Shapiro filters) on discrete $L_k$, as unguided smoothing corrupts $L''(z^*)$ and attenuates genuine fold geometry. Instead, fit discrete fluxes to shape-preserving splines or regularized local polynomials $L_h(z) \in C^2$ whose smoothing length scale is tied explicitly to grid resolution $\Delta z$.
* **Resolution-Aware Curvature Limiters:** Apply diagnostic limiters $\mathcal{L}(Ri_g, \partial_z Ri_g, \partial_{zz} Ri_g)$ strictly when a dimensionless curvature-resolution metric indicates under-resolution:

$$\mathcal{R}_{\text{fold}} = \frac{\frac{1}{2} \vert{}Ri_{g,zz}\vert{} (\Delta z)^2}{\vert{}Ri_g(z^*) - Ri_{\text{crit}}\vert{}}$$

Limiting activates only when $\mathcal{R}_{\text{fold}} \ge 1$, preserving fully resolved coordinate folds while preventing under-resolved grid cells from crossing $Ri_{\text{crit}}$.

---

**MEMORANDUM**

**TO:** Dick McNider & UAH Boundary Layer Meteorology Group

**SUBJECT:** Refined Numerical Formulations: Closure-Discretization Feedback and $C^2$ Reconstruction

---

**Executive Summary**

This revised memorandum updates our theoretical framework to distinguish physical coordinate folds from discretization-induced numerical feedback loops in planetary boundary layer (PBL) parameterizations.

---

**Key Formulations & Physical Corrections**

1. **Quadratic Sampling, Not Discontinuities:** A coordinate fold $\zeta_z(z^*) = 0$ generates a smooth quadratic variation $Ri_g(z^* \pm \Delta z) = Ri_g(z^*) + \frac{1}{2} Ri_{g,zz}(z^*) (\Delta z)^2$, not an artificial jump.
2. **Loop Gain Stability Criterion ($\mathcal{G}$):** Spurious quenching is defined as a closure–discretization instability driven by under-resolved curvature interacting with steep stability functions $\vert{}\mathcal{K}_h'\vert{} \gg 1$. The system runaway threshold occurs when the discrete feedback gain satisfies:

$$\mathcal{G} \sim \left(\frac{\partial Ri_g}{\partial \theta_z}\right) \left(\frac{\partial \theta_z}{\partial H}\right) \left(\frac{\partial H}{\partial K_h}\right) \left(\frac{\partial K_h}{\partial Ri_g}\right) \ge 1$$

1. **Diagnostic Taxonomy:** Separates legitimate coordinate geometry ($\zeta_z = 0, \zeta_{zz} \neq 0$) from under-resolved grid threshold crossing and gradient-estimation noise.

---

**Recommendations for Model Development**

* **Preserve Physical-Space Fluxes:** Retain governing transport equations in physical space $z$ rather than transforming differencing stencils to $\zeta$-space, avoiding local singularities near $\zeta_z \to 0$.
* **Shape-Preserving Spline Reconstructions:** Replace discrete grid filtering of $L_k$ with $C^2$-consistent functional fits $L_h(z)$ to evaluate derivative terms ($L', L''$) without distorting physical fold geometry.
* **Targeted Curvature Limiters:** Trigger diagnostic limiters using a curvature-resolution indicator $\mathcal{R}_{\text{fold}}$, ensuring fully resolved physical folds remain untouched.
