### Architectural Redesign: Physical, Constitutive, and Observational Separation

The revised GSPT Engine (`GSPTPhase2.jl`) separates profile processing into three explicit objects: **Coordinate Geometry ($G_\zeta$)**, **Constitutive Geometry ($G_{Ri}$)**, and **Observed Diagnostic ($O_{Ri}$)**.

```
                                 ProfileData Input
                        (z, u, v, θ_v, w'θ_v', u'w', v'w', e, K_m)
                                       │
                                       ▼
                       ┌───────────────────────────────┐
                       │   Coordinate Geometry G_ζ     │
                       │     (ζ,  ζ_z,  ζ_zz)          │
                       └───────────────┬───────────────┘
                                       │
                ┌──────────────────────┴──────────────────────┐
                ▼                                             ▼
┌───────────────────────────────┐             ┌───────────────────────────────┐
│   Constitutive Geometry G_Ri  │             │   Observation Diagnostic O_Ri │
│ (Ri_gspt, Ri_z, Ri_zz,        │             │ (Ri_obs, Ri_zz_obs,           │
│  C_const, C_coord, R_coord)   │             │  τ_Δz truncation, Δ_closure)  │
└───────────────┬───────────────┘             └───────────────┬───────────────┘
                │                                             │
                └──────────────────────┬──────────────────────┘
                                       ▼
                       ┌───────────────────────────────┐
                       │      Domain Comparison        │
                       │ Δ_geom,  Δ_closure,  Metrics  │
                       └───────────────────────────────┘

```

---

### Failure Mode Taxonomy & Model Improvement Criteria

By comparing direct coordinate biases $\Delta_{\text{geom}} = \mathcal{R}_{\text{coord}}^{\text{model}} - \mathcal{R}_{\text{coord}}^{\text{obs}}$ against closure residuals $\Delta_{\text{closure}} = Ri_{zz}^{\text{obs}} - Ri_{zz}^{\text{gspt}}$, modeling teams can isolate three distinct failure modes:

| Failure Mode | Geometric Condition ($\Delta_{\text{geom}}$) | Closure Residual ($\Delta_{\text{closure}}$) | Physical Mechanism | Corrective Modeling Action |
| --- | --- | --- | --- | --- |
| **I. Constitutive Failure** | $\Delta_{\text{geom}} \approx 0$ | $\vert{}\Delta_{\text{closure}}\vert{} \gg 0$ | Stability functions ($\phi_m, \phi_h$) fail to capture true local fluid shear thermodynamics. | Retune empirical stability parameters ($\beta_m, \beta_h$) or modify mixing length scales. |
| **II. Coordinate Distortion** | $\vert{}\Delta_{\text{geom}}\vert{} \gg 0$ | $\Delta_{\text{closure}} \approx 0$ | Model flux divergence ($L'(z)$) is too rigid, forcing incorrect coordinate stretching/compression. | Do **not** retune $\phi_m/\phi_h$. Modify non-local flux profiles or adjust vertical grid spacing near the jet. |
| **III. Coupled Breakdown** | $\vert{}\Delta_{\text{geom}}\vert{} \gg 0$ | $\vert{}\Delta_{\text{closure}}\vert{} \gg 0$ | Grid resolution is insufficient to resolve fold singularities alongside thermodynamic extinction. | Perform dynamic vertical grid adaptation ($\Delta z$) using $\vert{}\zeta_{zz}(z)\vert{}$ as the refinement indicator. |

---

### Methodological Refinements Applied

* **Independent Triple-Point ($C_{\text{TP}}$):** Replaced the midpoint proxy with an independent diffusivity extinction height ($z_K$, where $K_m \to K_{\min}$), ensuring $z_e, z_{ez}, z_K$ evaluate three distinct physical fields. Added `tke_fraction` parameterization for sensitivity testing.
* **Strict Mask Alignment:** Evaluated joint boolean vectors (`mask_scm`, `mask_les`) to eliminate silent `NaN` contamination during vertical bias and summary metric calculation ($\text{RMSE}$, $\text{MAE}$, $\text{median}\vert{}\Delta R\vert{}$).
* **Virtual Potential Temperature Scaling:** Standardized all Monin–Obukhov scale calculations to use virtual potential temperature ($\theta_v$) and kinematic virtual heat flux ($w'\theta_v'$).
* **Explicit Truncation Profiling ($\tau_{\Delta z}$):** Quantified Taylor discretization error $\tau_{\Delta z} = D_2 Ri_{\text{gspt}} - (C_{\text{const}} + C_{\text{coord}})$ to isolate grid-spacing sensitivity on irregular towers.