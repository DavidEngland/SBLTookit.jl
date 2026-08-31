# Coordinate Curvature Under Nieuwstadt Local Scaling: A Formal Treatment

## 1. Setup and Definitions

Let $z \in I = (0, z_{\text{top}})$ denote height above the surface within the stable boundary layer (SBL). Let $\tau(z) > 0$ denote the local vertical shear stress (momentum flux) and $H(z) \equiv \overline{w'\theta_v'}(z)$ denote the local kinematic virtual potential temperature flux, with $H(z) < 0$ under stable stratification. Let $\kappa$ be the von Kármán constant and $g/\theta_0$ the buoyancy parameter.

**Definition 1 (Local Obukhov length and inverse profile).**


$$L(z) \equiv \frac{-\tau(z)^{3/2}}{\kappa \dfrac{g}{\theta_0} H(z)}, \qquad \chi(z) \equiv \frac{1}{L(z)} = -\kappa \frac{g}{\theta_0} \frac{H(z)}{\tau(z)^{3/2}}$$

**Definition 2 (Local similarity coordinate).**


$$\zeta(z) = \frac{z}{L(z)} = z\,\chi(z)$$

**Standing Assumption (A1).** On any subinterval of interest, $\tau(z) > 0$ (turbulence is locally active), $H(z) < 0$, and $\chi \in C^2$, except possibly at isolated points explicitly treated in Lemma 4.

---

## 2. Coordinate Derivatives and Kinematics

**Proposition 1 (Coordinate derivatives in inverse-profile form).**
Under (A1), for all $z \in I$:


$$\zeta_z(z) = \chi(z) + z\,\chi'(z) \tag{1}$$

$$\zeta_{zz}(z) = 2\chi'(z) + z\,\chi''(z) \tag{2}$$

*Proof.* Immediate from Definition 2 by repeated application of the product rule:


$$\zeta_z = \frac{d}{dz}\big(z\chi\big) = \chi + z\chi'$$

$$\zeta_{zz} = \frac{d}{dz}\big(\chi + z\chi'\big) = \chi' + \big(\chi' + z\chi''\big) = 2\chi' + z\chi'' \quad \blacksquare$$

**Remark ( singularity removal).** Equations (1)–(2) are algebraically equivalent to the $L(z)$-based forms $\zeta_z = (1 - \zeta L')/L$ and $\zeta_{zz} = -(L''/L)\zeta - 2(L'/L)\zeta_z$. The $\chi$-parametrization removes the singularity as $L \to \infty$, making it well-suited for asymptotic analysis near neutral or collapse conditions.

Expanding $\chi'(z)/\chi(z) = (H'/H) - \frac{3}{2}(\tau'/\tau)$ allows $\zeta_z$ to be expressed directly in terms of relative flux gradients:


$$\zeta_z(z) = \chi(z) \left[ 1 + z \left( \frac{H'(z)}{H(z)} - \frac{3}{2}\frac{\tau'(z)}{\tau(z)} \right) \right] \tag{3}$$

---

## 3. Exclusion of Coordinate Folds and Differential Flux Divergence

Define a **coordinate fold point** $z_* \in I$ as a point where $\zeta_z(z_*) = 0$. It is **nondegenerate** if additionally $\zeta_{zz}(z_*) \neq 0$.

**Proposition 2 (Structural exclusion of folds under constant flux).**
Let $J \subset I$ be a subinterval on which $\chi$ is constant ($\chi \equiv \chi_0 \neq 0$, corresponding to a classical constant-flux surface layer). Then $\zeta_z(z) \neq 0$ for every $z \in J$; hence $J$ contains no coordinate fold points.

*Proof.* On $J$, $\chi' \equiv 0$, so by (1), $\zeta_z(z) = \chi_0$ for all $z \in J$. Since $\chi_0 \neq 0$, $\zeta_z$ never vanishes on $J$. $\blacksquare$

**Corollary 2.1 (Fold kinematics and differential divergence).**
Let $z_* > 0$ be a coordinate fold point under stable conditions ($\chi(z_*) > 0$). Then:


$$\chi'(z_*) = -\frac{\chi(z_*)}{z_*} < 0 \tag{4}$$


Equivalently, the fold condition is governed by differential flux divergence:


$$1 + z_* \left( \frac{H'(z_*)}{H(z_*)} - \frac{3}{2}\frac{\tau'(z_*)}{\tau(z_*)} \right) = 0 \tag{5}$$

*Proof.* Setting $\zeta_z(z_*) = 0$ in (1) yields $\chi(z_*) + z_*\chi'(z_*) = 0$, proving (4). Under assumption (A1), $\chi(z_*) > 0$ and $z_* > 0$, forcing $\chi'(z_*) < 0$. Substituting the expanded logarithmic derivative into $\zeta_z(z_*) = 0$ yields (5). $\blacksquare$

A coordinate fold in a turbulent, stable layer requires locally decreasing $\chi = 1/L$, driven by differential divergence between heat and momentum fluxes rather than thermal flux divergence alone.

---

## 4. Curvature Decomposition and the Fold Signature

Let $R(\cdot)$ be a smooth similarity function such that $Ri_g(z) = R(\zeta(z))$ with $R \in C^2$. High-curvature physical features ($\vert{}Ri_{g,zz}\vert{} \gg 0$) are distinguished from formal inflection points ($Ri_{g,zz} = 0$).

**Proposition 3 (Curvature decomposition).**


$$\frac{d Ri_g}{dz} = R'(\zeta)\,\zeta_z \tag{6}$$

$$\frac{d^2 Ri_g}{dz^2} = \underbrace{R''(\zeta)\,\zeta_z^{\,2}}_{\text{intrinsic similarity curvature}} + \underbrace{R'(\zeta)\,\zeta_{zz}}_{\text{coordinate-induced curvature}} \tag{7}$$

*Proof.* Direct application of the chain rule twice to $Ri_g(z) = R(\zeta(z))$. $\blacksquare$

**Theorem (Curvature signature of a nondegenerate fold).**
Let $Ri_g(z) = R(\zeta(z))$ with $R \in C^2$ and $\zeta \in C^2$. If $z_*$ is a nondegenerate coordinate fold ($\zeta_z(z_*) = 0, \zeta_{zz}(z_*) \neq 0$) and $R'(\zeta(z_*)) \neq 0$, then:


$$Ri_{g,z}(z_*) = 0 \qquad \text{and} \qquad Ri_{g,zz}(z_*) = R'(\zeta(z_*))\zeta_{zz}(z_*) \neq 0 \tag{8}$$

*Proof.* Evaluating (6) and (7) at $z_*$ where $\zeta_z(z_*) = 0$ eliminates the intrinsic similarity term ($R''\zeta_z^2 = 0$), leaving $Ri_{g,zz}(z_*) = R'(\zeta(z_*))\zeta_{zz}(z_*)$. Because $R'(\zeta(z_*)) \neq 0$ and $\zeta_{zz}(z_*) \neq 0$, their product is non-zero. $\blacksquare$

A nondegenerate coordinate fold guarantees a stationary point with nonzero physical-space curvature in $Ri_g(z)$, irrespective of the intrinsic curvature $R''(\zeta)$.

**Corollary 3.1 (Quantitative curvature decomposition ratio).**
To quantify whether an observed physical-space curvature feature ("knee", $\vert{}Ri_{g,zz}\vert{} \gg 0$) stems from similarity closure geometry or coordinate mapping, define:


$$\mathcal{C}_{\text{coord}}(z) \equiv \frac{\vert{}R'(\zeta)\,\zeta_{zz}\vert{}}{\vert{}R''(\zeta)\,\zeta_z^{\,2}\vert{} + \vert{}R'(\zeta)\,\zeta_{zz}\vert{} + \varepsilon_C} \qquad (\varepsilon_C > 0)$$

* $\mathcal{C}_{\text{coord}} \approx 1$: Coordinate-induced curvature dominates; apparent physical features are driven by non-linear coordinate mapping ($\zeta_{zz}$).
* $\mathcal{C}_{\text{coord}} \approx 0$: Intrinsic similarity curvature dominates; physical features reflect structural curvature in the scaling function $R(\zeta)$.
* Intermediate values: Both mechanisms contribute comparably.

---

## 5. Regularity at Turbulence-Collapse Heights

Near a turbulence-collapse height $z_c$, where $\tau(z) \to 0$ and $H(z) \to 0$ as $z \to z_c^-$, $\chi(z)$ forms a $0/0$ indeterminate limit. Assume smooth, twice-differentiable asymptotic expansions:


$$\tau(z) = A(z_c - z)^a [1 + o(1)], \qquad H(z) = -B(z_c - z)^b [1 + o(1)] \qquad (A, B, a, b > 0)$$


The local scaling parameter behaves as $\chi(z) \sim C(z_c - z)^p$ where $C = \kappa \frac{g}{\theta_0} \frac{B}{A^{3/2}}$ and $p \equiv b - \frac{3}{2}a$.

**Lemma 4 (Asymptotic regularity spectrum of $\chi$ at $z_c$).**
Under controlled asymptotic differentiation, the regularity of $\chi(z)$ as $z \to z_c^-$ is governed by $p$:

| Exponent $p = b - \frac{3}{2}a$ | Regularity and Limiting Behavior as $z \to z_c^-$ |
| --- | --- |
| $p < 0$ | Singular coordinate mapping: $\chi \to \infty$ |
| $p = 0$ | Discontinuous limit: $\chi \to C \neq 0$ |
| $0 < p < 1$ | $C^0$ regular: $\chi \to 0$, but $\chi'$ is singular |
| $1 \le p < 2$ | $C^1$ regular: $\chi' \to 0$, but $\chi''$ is singular |
| $p = 2$ | $C^2$ bounded: $\chi'' \to 2C \neq 0$ (finite non-zero limit) |
| $p > 2$ | $C^2$ vanishing: $\chi, \chi', \chi'' \to 0$ |

*Remark (Near-neutral conditions).* For pure near-neutral decay ($a = 0 \implies p = b$), $b \ge 2$ is required for $\zeta_{zz}$ to remain bounded at $z_c$, while $b > 2$ is required for $\chi'' \to 0$. Linear flux decay ($a = 1, b = 1 \implies p = -0.5$) breaks $C^2$ regularity, generating a coordinate singularity at $z_c$.

---

## 6. Summary of Formal Results

| Result | Mathematical Statement | Operational / Physical Interpretation |
| --- | --- | --- |
| **Prop. 1** | $\zeta_z = \chi + z\chi'$, $\zeta_{zz} = 2\chi' + z\chi''$ | Singularity-free derivative formulation via $\chi(z) = 1/L(z)$ |
| **Prop. 2 / Cor. 2.1** | $\zeta_z(z_*) = 0 \implies 1 + z_*\left(\frac{H'}{H} - \frac{3}{2}\frac{\tau'}{\tau}\right) = 0$ | Folds are excluded in constant-flux layers and require differential flux divergence |
| **Theorem** | $\zeta_z(z_*) = 0 \implies Ri_{g,zz}(z_*) = R'\zeta_{zz} \neq 0$ | Nondegenerate folds enforce physical curvature independently of $R''$ |
| **Cor. 3.1** | $\mathcal{C}_{\text{coord}} = \frac{\Vert{}R'\zeta_{zz}\Vert{}}{\Vert{}R''\zeta_z^2\Vert{} + \Vert{}R'\zeta_{zz}\Vert{} + \varepsilon_C}$ | Diagnostic metric separating coordinate mapping from intrinsic closure curvature |
| **Lemma 4** | $p = b - \frac{3}{2}a \ge 2$ | Threshold for bounded $C^2$ coordinate extensions near collapse heights |