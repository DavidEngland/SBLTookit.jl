# Geometric Analysis of Richardson Profile Curvature and Coordinate Folds in the Stable Boundary Layer

Under local scaling theory (Nieuwstadt 1984), physical height $z \in \mathbb{R}^+$ maps to dimensionless local similarity space via $\zeta(z) = z / L(z)$, where $L(z) > 0$ varies continuously across the nocturnal stable boundary layer (SBL). The observed gradient Richardson number profile $Ri_g(z) = R(\zeta(z))$ decomposes its total spatial curvature into exact intrinsic thermodynamic and coordinate-mapping components:

$$Ri_{g,zz}(z) = \underbrace{R''(\zeta)\zeta_z^2}_{C_{\text{closure}}(z)} + \underbrace{R'(\zeta)\zeta_{zz}}_{C_{\text{mapping}}(z)}$$

where spatial coordinate derivatives are given by:

$$\zeta_z = \frac{L(z) - z L'(z)}{L(z)^2}, \qquad \zeta_{zz} = \frac{-z L''(z) L(z) - 2 L'(z)\left(L(z) - z L'(z)\right)}{L(z)^3}$$

---

**Spatial Coordinate Fold Taxonomy**

The condition $\zeta_z(z^*) = 0$ defines a critical point of the composite mapping $\zeta(z) = z/L(z)$—requiring $L'(z^*) = L(z^*)/z^*$—rather than an extremum of $L(z)$ itself ($L'(z^*) = 0$). Spatial fold structures are classified into three distinct hierarchical tiers:

* **Coordinate Critical Point ($\zeta_z(z^*) = 0$):** A height $z^*$ where $L(z^*) - z^* L'(z^*) = 0$. The closure curvature vanishes identically ($C_{\text{closure}}(z^*) = 0$).
* **Nondegenerate Coordinate Fold ($\zeta_z(z^*) = 0, \, \zeta_{zz}(z^*) \neq 0$):** Requires $z^* L''(z^*) \neq 0$ (assuming $L(z^*) \neq 0$), establishing a non-zero geometric turning structure in the mapping space.
* **Visible Richardson-Curvature Fold ($\zeta_z(z^*) = 0, \, \zeta_{zz}(z^*) \neq 0, \, R'(\zeta^*) \neq 0$):** Requires non-zero local closure sensitivity, ensuring the coordinate fold projects cleanly into observable physical profile curvature:

$$\left. Ri_{g,zz} \right\vert{}_{z^*} = -R'(\zeta^*) \frac{z^* L''(z^*)}{L(z^*)^2} \neq 0$$

---

**Businger–Dyer Closure Evaluation ($\beta_m = \beta_h = \beta > 0$)**

Assuming standard positive stability coefficients ($\beta > 0$), the closure $R(\zeta) = \frac{\zeta(Pr_0^{-1} + \beta\zeta)}{(1 + \beta\zeta)^2}$ yields:

$$R'(\zeta) = \frac{Pr_0^{-1} + \beta(2 - Pr_0^{-1})\zeta}{(1 + \beta\zeta)^3}, \qquad R''(\zeta) = \frac{2\beta\left[(Pr_0 - 2) - \beta\zeta(2 Pr_0 - 1)\right]}{Pr_0(1 + \beta\zeta)^4}$$

* **Canonical Regime ($0.5 \le Pr_0 \le 2$):** $R''(\zeta) < 0$ strictly for all $\zeta > 0$. The closure contains no intrinsic positive-$\zeta$ curvature reversals. Any observed physical inflection ($Ri_{g,zz} = 0$) away from a fold is co-produced by the exact spatial balance $C_{\text{closure}} + C_{\text{mapping}} = 0$.
* **Neutral Baseline ($Pr_0 = 1$):** Reduces to the Webb profile $R(\zeta) = \frac{\zeta}{1 + \beta\zeta}$, yielding a visible fold curvature of:

$$\left. Ri_{g,zz} \right\vert{}_{z^*} = -\frac{z^* L''(z^*)}{L(z^*)^2 (1 + \beta\zeta^*)^2}$$

* **Anomalous Regimes ($Pr_0 > 2$ or $Pr_0 < 0.5$):** Admits an intrinsic similarity inflection point at $\zeta_{\text{inf}} = \frac{Pr_0 - 2}{\beta(2 Pr_0 - 1)}$.

---

**Implications for SBL Diagnostics and Numerical Models**

* **Physical Interpretation:** Profile knees in observational tower or remote-sensing data reflect curvature projection by a folded similarity coordinate rather than intrinsic thermodynamic singularities of $R(\zeta)$. At $z^*$, curvature is purely mapping-controlled; away from $z^*$, coordinate geometry co-shapes the observed structure alongside thermodynamic closure sensitivity.
* **Model Regularization:** Evaluating $C_{\text{mapping}}$ requires $L''(z)$, a derivative of a flux-derived field highly sensitive to noise. Boundary layer schemes and diagnostic algorithms must utilize smooth functional fits or regularized derivative operators to prevent false turbulent quenching near shear zones.

---

**MEMORANDUM**

**TO:** Dick McNider & UAH Boundary Layer Meteorology Group

**SUBJECT:** Updated Framework: Profile Curvature Decomposition and Fold Taxonomy in the SBL

---

**Executive Summary**

This updated draft presents a continuous differential-geometric framework for isolating intrinsic turbulence-state transitions from coordinate-induced mapping geometry in observed gradient Richardson number profiles ($Ri_g(z)$).

Using Nieuwstadt (1984) local scaling ($\zeta(z) = z / L(z)$), spatial profile curvature decomposes into:

$$Ri_{g,zz}(z) = \underbrace{R''(\zeta)\zeta_z^2}_{C_{\text{closure}}(z)} + \underbrace{R'(\zeta)\zeta_{zz}}_{C_{\text{mapping}}(z)}$$

---

**Core Refinements & Taxonomy**

* **Mapping Criticality vs. Profile Extrema:** The fold condition $\zeta_z(z^*) = 0 \iff L'(z^*) = L(z^*)/z^*$ defines a critical point of the composite coordinate transformation $z \mapsto \zeta$, rather than a local extremum of $L(z)$ itself ($L'(z^*) = 0$).
* **Three-Tier Nondegeneracy Classification:**

1. *Coordinate Critical Point:* $\zeta_z(z^*) = 0 \implies C_{\text{closure}}(z^*) = 0$.
2. *Nondegenerate Coordinate Fold:* $\zeta_z(z^*) = 0$ and $L''(z^*) \neq 0 \implies \zeta_{zz}(z^*) \neq 0$.
3. *Visible Richardson-Curvature Fold:* $\zeta_z(z^*) = 0$, $L''(z^*) \neq 0$, and $R'(\zeta^*) \neq 0 \implies Ri_{g,zz}(z^*) = -R'(\zeta^*)\frac{z^* L''(z^*)}{L(z^*)^2} \neq 0$.

* **Closure Non-Reversal ($\beta > 0$):** In the canonical regime ($0.5 \le Pr_0 \le 2$), $R''(\zeta) < 0$ strictly for all $\zeta > 0$. Profile knees away from $z^*$ represent a superposition balance ($C_{\text{closure}} + C_{\text{mapping}} = 0$), whereas curvature strictly at $z^*$ is purely mapping-generated.

---

**Diagnostic Recommendations for UAH Group**

1. **Observational Diagnostics:** Field profile inflections near Low-Level Jets (LLJs) or inversion tops should be evaluated as curvature projections by a folded similarity coordinate rather than assumed turbulence breakdowns.
2. **Noise Regularization:** Because $C_{\text{mapping}}$ depends on $L''(z)$, observational and numerical workflows must apply smooth functional representations of $L(z)$ to prevent noise amplification when calculating second derivatives.
