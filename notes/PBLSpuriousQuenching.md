Spurious quenching in planetary boundary layer (PBL) parameterizations occurs when under-resolved similarity profile curvature interacts with steep non-linear stability functions, initiating an unphysical positive feedback loop that collapses turbulent heat transport.

## Theoretical Background

* **Pathology Mechanism:** In stable boundary layers, spatial variation around coordinate folds ($\zeta_z \to 0$) generates local quadratic curvature in the Richardson profile: $Ri_g(z^* \pm \Delta z) \approx Ri_g(z^*) + \frac{1}{2} Ri_{g,zz}(z^*) (\Delta z)^2$. When vertical grid spacing $\Delta z$ cannot resolve this curvature, discrete $Ri_g$ artificially crosses the critical threshold $Ri_{\text{crit}} = 1/\beta$.
* **Feedback Loop Gain ($\mathcal{G}$):** A discrete threshold crossing collapses eddy diffusivity $K_h$, halting downward sensible heat flux $H = -\rho c_p K_h \theta_z$. This isolates the surface, rapidly cooling the lowest layer and steepening potential temperature gradients $\theta_z$. The system runaway instability is governed by loop gain:

$$\mathcal{G} = \underbrace{\left(\frac{\partial Ri_g}{\partial \theta_z}\right) \left(\frac{\partial \theta_z}{\partial H}\right) \left(\frac{\partial H}{\partial K_h}\right) \mathcal{K}_h'(Ri_g^{\text{reg}})}_{\mathcal{G}_{\text{physical}}} \cdot \mathcal{D}_\alpha(Ri_g)$$

* **Diagnostic Regularization:** The algorithm regularizes $Ri_g$ using a $C^1$-continuous algebraic mapping before evaluating $K_h = \mathcal{K}_h(Ri_g^{\text{reg}})$. This inserts a damping factor $\mathcal{D}_\alpha(Ri_g) \in \left[\frac{1}{1+\alpha}, \, 1\right)$ that strictly bounds sensitivity at threshold without modifying prognostic physical-space transport equations.

## Algorithm Parameter Reference

| Parameter | Symbol | Standard Value | Role & Dimensional Constraints |
| --- | --- | --- | --- |
| BD Constant | $\beta$ | $5.0$ | Dimensionless; sets stability decay slope |
| Critical Richardson | $Ri_{\text{crit}}$ | $0.20$ | Dimensionless; asymptotic bound $1/\beta$ |
| Softening Term | $\epsilon$ | $0.05$ | Dimensionless; guarantees $C > 0$ globally |
| Damping Control | $\alpha$ | $2.0$ | Dimensionless; bounds sensitivity $1/(1+\alpha) \le \mathcal{D}_\alpha < 1$ |
| Spline Penalty | $\gamma$ | $1.0$ | Dimensionless; scales penalty parameter $\lambda = \gamma \, \Delta z$ ($\text{m}^1$) |

## Step-by-Step Implementation Pipeline

1. **$C^2$ Obukhov Profile Reconstruction:** At each column timestep, fit discrete Obukhov lengths $L_k$ to a shape-preserving cubic smoothing spline $L_h(z)$ minimizing:

$$\min_{L_h} \sum_{k=1}^N \left(\frac{L_h(z_k) - L_k}{\sigma_k}\right)^2 + \lambda \int_0^h \left[L_h''(z)\right]^2 \, \mathrm{d}z, \qquad \text{where } \lambda = \gamma \, \Delta z$$

1. **Analytic Differential Mapping:** Evaluate smooth derivatives $L_h'(z_k)$ and $L_h''(z_k)$ to compute continuous geometric derivatives at cell center $z_k$:

$$\zeta = \frac{z_k}{L_h}, \qquad \zeta_z = \frac{L_h - z_k L_h'}{L_h^2}, \qquad \zeta_{zz} = \frac{-z_k L_h'' L_h - 2 L_h'(L_h - z_k L_h')}{L_h^3}$$

1. **Profile Curvature & Resolution Scale:** Evaluate total physical profile curvature $Ri_{g,zz}$ using Businger–Dyer derivatives:

$$Ri_{g,zz} = \frac{-2\beta}{(1 + \beta \zeta)^3} \, \zeta_z^2 + \frac{1}{(1 + \beta \zeta)^2} \, \zeta_{zz}$$

$$C = \frac{1}{2} \vert{}Ri_{g,zz}\vert{} (\Delta z)^2 + \epsilon \, Ri_{\text{crit}}$$

1. **Smooth Diagnostic Regularization:** Apply the $C^1$ smooth algebraic limiter to map raw diagnostic $Ri_g$ to regularized $Ri_g^{\text{reg}}$:

$$Ri_g^{\text{reg}} = Ri_{\text{crit}} + (Ri_g - Ri_{\text{crit}}) \left[ \frac{\vert{}Ri_g - Ri_{\text{crit}}\vert{} + C}{\vert{}Ri_g - Ri_{\text{crit}}\vert{} + (1 + \alpha) C} \right]$$

1. **Eddy Diffusivity & Physical Flux Step:** Update mixing diffusivities and physical-space turbulent heat transport:

$$K_h = \mathcal{K}_h(Ri_g^{\text{reg}}), \qquad H = -\rho c_p K_h \frac{\partial \theta}{\partial z}$$
