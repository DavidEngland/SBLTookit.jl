The **Inflection Masking Illusion** occurs when a vertical gradient Richardson number ($Ri_g$) sounding appears locally linear ($Ri_{zz} \approx 0$). Generalized Similarity Profile Theory (GSPT) proves that a straight physical profile is a geometric artifact where physics-driven thermodynamic cooling ($C_{\text{const}}$) is actively neutralized by non-uniform coordinate-stretching geometry ($C_{\text{coord}}$):

$$\frac{d^2 Ri_g}{dz^2} = \underbrace{R''(\zeta)\zeta_z^2}_{C_{\text{const}}} + \underbrace{R'(\zeta)\zeta_{zz}}_{C_{\text{coord}}} + \mathcal{E}_{\Delta z}$$

**Quantitative Partitioning at $z \approx 10\text{ m}$ (Hours 6–12)**

* **Intrinsic Thermodynamic Curvature ($C_{\text{const}}$):** $-0.0040\text{ m}^{-2}$ (smooth stabilization of lower surface layers)
* **Coordinate Mapping Curvature ($C_{\text{coord}}$):** $+0.0036\text{ m}^{-2}$ (flux-divergent grid deformation)
* **Observed Physical Curvature ($Ri_{zz}$):** $-0.0004\text{ m}^{-2}$ (near-zero net curvature due to offsetting terms)

**Altitudinal Regime Comparison**

| Boundary Layer Region | Stability Coordinate | Curvature Partition | Physical & Geometric Dynamics |
| --- | --- | --- | --- |
| **Near-Surface Layer** | $\zeta \to 0$ | Term Separation ($C_{\text{const}} \ll C_{\text{coord}}$) | Thermodynamic non-linearity operates at full strength without saturation dampening ($R''(\zeta) \approx -2\beta$). Intrinsic curvature dominates, exposing true physical profile bending. |
| **Upper SBL / Jet Altitude** | $\zeta \gg 0.5$ | Perfect Cancellation ($C_{\text{coord}} \approx -C_{\text{const}}$) | Profile flattens near its saturated ceiling ($Ri_{\text{sat}} = 1/\beta \approx 0.20$). Cubic saturation $(1 + \beta \zeta)^{-3}$ dampens both terms, maintaining the straight profile illusion. |

**NWP & Parameterization Hazards**

* **Compensating Errors:** Single-Column Models (SCMs) process physical-space gradients directly and cannot distinguish between neutral stratification and active inflection masking under strong flux divergence.
* **Degraded Mixing Physics:** Retuning empirical Monin–Obukhov stability constants ($\beta_m, \beta_h$) to force constitutive equations to match geometrically masked profiles introduces severe over-dampening in un-warped, homogeneous regimes.
