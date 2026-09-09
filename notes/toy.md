Throughout, $\zeta(z) \equiv z/L(z)$ is the Monin–Obukhov stability parameter, and $\phi_h, \phi_m$
are the similarity functions for heat and momentum. Analysis restricted to the stable branch
$\zeta \ge 0$ unless noted.

* **Richardson Mapping Function $R(\zeta)$:**

$$R(\zeta) = \zeta \frac{\phi_h(\zeta)}{\phi_m^2(\zeta)}, \qquad Ri_g(z) = R(\zeta(z))$$

* **Height-Resolved Curvature Decomposition:**

$$\frac{d^2 Ri_g}{dz^2} = \underbrace{R''(\zeta)\,\zeta_z^2}_{C_{\phi}(z)} + \underbrace{R'(\zeta)\,\zeta_{zz}}_{C_{\text{coord}}(z)} + \mathcal{E}_{\Delta z}$$

  where $\mathcal{E}_{\Delta z}$ is the finite-difference truncation error incurred when $d^2Ri_g/dz^2$
  is estimated numerically on a grid of spacing $\Delta z$ (it vanishes for the exact analytic
  derivative). $C_\phi$ captures curvature intrinsic to the stability function's shape;
  $C_{\text{coord}}$ captures curvature injected purely by the nonlinearity of the $z\to\zeta$
  coordinate map.

* **Coordinate Jacobian $\zeta_z$:**

$$\zeta_z = \frac{L(z) - z L'(z)}{L(z)^2} = \frac{1 - \zeta L'(z)}{L(z)}$$

* **Linear Obukhov Test (No-Fold Gate):**
If $L(z) = L_0 + az$ with $L_0 > 0$ and $a \ge 0$ (so $L(z) \neq 0$ for all $z \ge 0$):

$$\zeta_z = \frac{L_0}{(L_0 + az)^2} > 0 \quad \forall z \ge 0$$

  No finite coordinate turning points can exist. (If $a<0$, this holds only up to
  $z < -L_0/a$, beyond which $L\to0$ and $\zeta$ itself becomes singular.)

* **Nondegenerate Fold Condition:**
At a fold locus $z^* > 0$ where $L(z^*) \neq 0$:

$$\zeta_z(z^*) = 0 \iff L(z^*) = z^* L'(z^*), \quad L''(z^*) \neq 0$$

$$\zeta_{zz}(z^*) = -\frac{z^* L''(z^*)}{L(z^*)^2}$$

* **Curvature at a Fold (exact, up to discretization error):**

$$\frac{d^2Ri_g}{dz^2}\bigg|_{z=z^*} = R'(\zeta^*)\,\zeta_{zz}(z^*) + \mathcal{E}_{\Delta z}$$

  Strictly coordinate-driven at the fold, since $C_\phi(z^*) = 0$ identically.

---

Applying the rational, saturated stability function $R(\zeta) = \dfrac{\zeta}{1 + \beta \zeta}$
($\beta > 0$: high-stability saturation parameter) to $Ri_g(z) = R(\zeta(z))$ yields closed-form
curvature and fold expressions.

**1. Derivative Operators of $R(\zeta)$**

$$R'(\zeta) = \frac{1}{(1 + \beta \zeta)^2}, \qquad R''(\zeta) = -\frac{2\beta}{(1 + \beta \zeta)^3}$$

**2. Curvature Decomposition**

$$C_\phi(z) = -\frac{2\beta}{\left(1 + \beta \zeta(z)\right)^3} \left[ \frac{1 - \zeta(z) L'(z)}{L(z)} \right]^2$$

$$C_{\text{coord}}(z) = \frac{1}{\left(1 + \beta \zeta(z)\right)^2} \left[ -\frac{\zeta(z) L''(z)}{L(z)} - \frac{2 L'(z)}{L(z)} \zeta_z(z) \right]$$

* **Intrinsic Stabilization ($C_\phi$):** For $\zeta \ge 0$ and $\beta>0$, $1+\beta\zeta>0$, so
  $C_\phi(z) \le 0$ everywhere on the stable branch — smooth, monotonic thermodynamic
  stabilization with no sign changes. (This sign flips past the function's pole at
  $\zeta = -1/\beta$, which is outside the physical domain considered here.)
* **Saturation Dampening:** The $(1+\beta\zeta)^3$ denominator suppresses $C_\phi$ as
  $\zeta \to \infty$, preventing unphysical intrinsic curvature blowup under strongly
  stable conditions.

**3. Fold & Turning Point Analysis at $z = z^*$**

A spatial coordinate fold occurs where the coordinate Jacobian vanishes: $\zeta_z(z^*) = 0$.

* **Turning Condition:**

$$\zeta_z(z^*) = 0 \iff L(z^*) = z^* L'(z^*) \implies \zeta(z^*) = \frac{1}{L'(z^*)}$$

* **Nondegeneracy Condition:** the fold is non-degenerate iff $\zeta_{zz}(z^*) \neq 0$,
  which requires $L''(z^*) \neq 0$ (and $z^* \neq 0$).

* **Curvature Partition at $z^*$:**

$$C_\phi(z^*) = -\frac{2\beta}{(1 + \beta \zeta^*)^3}\,(0)^2 = 0$$

$$C_{\text{coord}}(z^*) = -\frac{z^* L''(z^*)}{L(z^*)^2 \left(1 + \beta \dfrac{z^*}{L(z^*)}\right)^2}$$

$$\frac{d^2 Ri_g}{dz^2}\bigg|_{z=z^*} = C_{\text{coord}}(z^*)$$

**Key Analytical Insight**

At any turning point $z^*$, $C_\phi$ vanishes identically because $\zeta_z(z^*) = 0$.
Consequently, **all of the observed profile curvature at a fold is coordinate-driven**
($C_{\text{coord}}$) — a "Fold Illusion": an apparent kink in $Ri_g(z)$ that reflects the
geometry of the $z \to \zeta$ map rather than any genuine change in atmospheric
stratification physics.

---

That cancellation in your simulations directly validates a central GSPT prediction: the **Inflection Masking Illusion**. When an observed Richardson profile appears locally linear ($Ri_{zz} \approx 0$), it is rarely due to trivial local physics, but rather coordinate-stretching geometry ($C_{\text{coord}} > 0$) actively neutralizing intrinsic thermodynamic curvature ($C_{\text{const}} < 0$).

**Mechanics of Upper-Boundary Cancellation vs. Surface Divergence**

* **Upper SBL Regime ($\zeta \gg 0.5$): Perfect Cancellation**
As stability deepens aloft, $Ri_g(z)$ asymptotically approaches its saturated ceiling ($Ri_{\text{sat}} = 1/\beta \approx 0.20$). Because the physical profile flattens ($Ri_{zz} \approx 0$), the curvature partition collapses to:

$$C_{\text{coord}} \approx -C_{\text{const}}$$

The cubic saturation factor $(1 + \beta \zeta)^{-3}$ dampens both terms, but their near-exact opposite signs create the visual illusion of a straight, un-warped profile at the top of the SBL.

* **Near-Surface Layer ($\zeta \to 0$): Term Separation**
Closer to the ground, $\zeta$ is small, so $R''(\zeta) \approx -2\beta$ operates at full strength without heavy saturation dampening. Intrinsic thermodynamic curvature ($C_{\text{const}}$) reaches its largest negative magnitude. Because $\zeta_z$ and $\zeta_{zz}$ scale differently with surface boundary conditions, $C_{\text{const}}$ and $C_{\text{coord}}$ separate, yielding a non-zero net $Ri_{zz}$ (a physically curved surface layer).
* **Why the Residual ($\mathcal{E}_{\Delta z}$) Remains Small**
The discretization residual $\mathcal{E}_{\Delta z}$ stays minimal across all levels because Track A pre-diagnostic log-height splines ($\xi = \ln(z/z_0)$) eliminate finite-difference stencil truncation errors. This confirms that the cancellation you are seeing is a true continuous feature of the governing mapping $R(\zeta(z))$, not a numerical artifact of the grid.
