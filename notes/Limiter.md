In Computational Fluid Dynamics (CFD) and Numerical Weather Prediction (NWP), limiters prevent spurious numerical oscillations around steep spatial gradients. While classical advective limiters (e.g., Total Variation Diminishing or TVD schemes) operate directly on prognostic scalar fluxes, this diagnostic limiter regularizes the gradient Richardson number ($Ri_g$) before it evaluates non-linear eddy diffusivity functions, suppressing runaway nocturnal cooling while leaving physical transport equations intact.

**Theoretical Background: Advective vs. Diagnostic Limiters**

* **Classical Advective Limiters:** Standard schemes (such as van Leer, Superbee, or minmod) constrain spatial advective fluxes $\mathbf{u} \cdot \nabla \phi$ to maintain monotonicity, preventing unphysical numerical overshoots near shocks or steep interfaces.
* **The Diagnostic Stability Bottleneck:** In planetary boundary layer (PBL) schemes, runaway cooling stems from the *diagnostic* pipeline rather than advection. When physical coordinate curvature $Ri_{g,zz}$ is under-resolved by grid spacing $\Delta z$, the raw discrete diagnostic $Ri_g$ jumps over $Ri_{\text{crit}} \approx 0.25$, dropping mixing diffusivity $K_h$ to a background floor and decoupling the surface.
* **The $C^\infty$ Algebraic Regularizer:** Rather than clipping physical fluxes after the fact, this diagnostic gate applies a smooth transformation to $Ri_g$ prior to evaluating $K_h = \mathcal{K}_h(Ri_g)$. This bounds the feedback loop gain $\vert{}\mathcal{G}\vert{} < 1$ without altering physical transport equations or introducing artificial derivative spikes.

---

**Step-by-Step Model Implementation Guide**

**Phase 1: Profile Reconstruction & Mapping Geometry**
At each column evaluation and time step:

1. Reconstruct a $C^2$-continuous Obukhov length field $L_h(z)$ from discrete grid values $L_k$ using a shape-preserving cubic smoothing spline with penalty parameter $\lambda = \gamma \, \Delta z$ (where $\gamma = 1.0$).
2. Evaluate continuous derivatives $L_h'(z_k)$ and $L_h''(z_k)$ analytically at each cell center $z_k$.
3. Compute mapping derivatives:

$$\zeta_z = \frac{L_h - z_k L_h'}{L_h^2}, \qquad \zeta_{zz} = \frac{-z_k L_h'' L_h - 2 L_h'(L_h - z_k L_h')}{L_h^3}$$

**Phase 2: Curvature Decomposition & Under-Resolution Metric**

1. Evaluate closure derivatives for local similarity coordinate $\zeta = z_k / L_h$. For Businger–Dyer ($Pr_0 = 1, \, \beta = 5.0$):

$$R'(\zeta) = \frac{1}{(1 + \beta \zeta)^2}, \qquad R''(\zeta) = \frac{-2 \beta}{(1 + \beta \zeta)^3}$$

1. Compute total physical profile curvature:

$$Ri_{g,zz} = R''(\zeta) \zeta_z^2 + R'(\zeta) \zeta_{zz}$$

1. Calculate the curvature resolution scale $C$:

$$C = \frac{1}{2} \vert{}Ri_{g,zz}\vert{} (\Delta z)^2 + \epsilon \, Ri_{\text{crit}}$$

**Phase 3: Smooth Regularization & Flux Update**

1. Map raw $Ri_g$ to regularized $Ri_g^{\text{reg}}$:

$$Ri_g^{\text{reg}} = Ri_{\text{crit}} + (Ri_g - Ri_{\text{crit}}) \left[ \frac{\vert{}Ri_g - Ri_{\text{crit}}\vert{} + C}{\vert{}Ri_g - Ri_{\text{crit}}\vert{} + (1 + \alpha) C} \right]$$

1. Pass $Ri_g^{\text{reg}}$ into the mixing scheme lookup: $K_h = \mathcal{K}_h(Ri_g^{\text{reg}})$.
2. Update vertical heat transport in physical space: $H = -\rho c_p K_h \frac{\partial \theta}{\partial z}$.

---

**Implementation Parameter Reference**

| Parameter | Description | Standard Value | Function in Algorithm |
| --- | --- | --- | --- |
| $\beta$ | BD Stability Constant | $5.0$ | Controls closure decay rate in $R'(\zeta)$ and $R''(\zeta)$ |
| $Ri_{\text{crit}}$ | Critical Richardson Number | $0.25$ | Threshold for turbulent diffusivity cutoff |
| $\epsilon$ | Denominator Softening | $0.05$ | Prevents metric zero-division near threshold |
| $\alpha$ | Relaxation Gain Control | $2.0$ | Enforces bounded sensitivity $\frac{\partial Ri_g^{\text{reg}}}{\partial Ri_g} \le \frac{1}{1 + \alpha}$ |
| $\gamma$ | Spline Penalty Factor | $1.0$ | Dimensions smoothing parameter $\lambda = \gamma \, \Delta z$ to grid spacing |
