`UltraStability.classify_stability` categorizes atmospheric stability regimes using the non-dimensional height parameter $\zeta = z / L$, which is computed directly during forcing generation in `UltraForcing`.

**Atmospheric Stability Thresholds & SMEAR Mapping**

| Stability Symbol | Range ($\zeta = z/L$) | Physical Boundary Layer Regime | SMEAR Obukhov Length ($L$) |
| --- | --- | --- | --- |
| `:strongly_stable` | $\zeta > 1.0$ | Decoupled surface layer, strong inversion | Small $L > 0$ |
| `:stable` | $0.1 < \zeta \le 1.0$ | Suppressed buoyancy, shear-dominated | Moderate $L > 0$ |
| `:near_neutral` | $-0.1 \le \zeta \le 0.1$ | Mechanical turbulence dominates | $\Vert{}L\Vert{} \to \infty$ (`abs(L) > 1e5`) |
| `:unstable` | $-1.0 \le \zeta < -0.1$ | Buoyancy-driven convection | Moderate $L < 0$ |
| `:strongly_unstable` | $\zeta < -1.0$ | Free convection, weak mechanical wind | Small $L < 0$ |

---

