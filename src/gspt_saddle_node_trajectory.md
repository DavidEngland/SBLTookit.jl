This technical reference details the slow-fast Stable Boundary Layer (SBL) dynamical framework for boundary-layer meteorologists, atmospheric physicists, and turbulence modelers.

---

**Governing Dynamical System**

The model describes the non-equilibrium interaction between turbulent kinetic energy ($e$) and vertical wind shear ($S$) under strong atmospheric stratification. The parameter $\epsilon = 0.05$ governs time-scale separation: TKE adjusts rapidly ($O(\epsilon)$ minutes) relative to environmental shear evolution ($O(1)$ hours).

$$\epsilon \frac{de}{dt} = F(e, S) = \underbrace{l_0 e S_{\text{eff}}^2}_{\text{Mechanical Production}} - \underbrace{B_{0,\max} \frac{e^2}{e^2 + \delta_{\text{reg}}^2}}_{\text{Buoyant Suppression}} - \underbrace{\frac{e^3}{l_0 \left(1 + \frac{\beta N^2}{S_{\text{eff}}^2}\right)}}_{\text{Dissipation}}$$

$$\frac{dS}{dt} = G(e, S) = \underbrace{G_0}_{\text{Geostrophic Forcing}} - \underbrace{\gamma_s e S}_{\text{Turbulent Mixing}} - \underbrace{r_s S}_{\text{Background Damping}}$$

where $S_{\text{eff}} = \sqrt{S^2 + \eta_S^2}$ provides a $C^\infty$ smooth formulation of shear near $S = 0$.

---

**Model Parameters & Meteorological Mapping**

| Parameter | Value / Units | Atmospheric & GSPT Meaning |
| --- | --- | --- |
| $e$ | $\mathrm{m^2\,s^{-2}}$ | Fast variable: Turbulent Kinetic Energy (TKE) intensity |
| $S$ | $\mathrm{s^{-1}}$ | Slow variable: Vertical wind shear magnitude ($\partial U / \partial z$) |
| $\epsilon$ | $0.05$ | Time-scale ratio ($\tau_{\text{TKE}} / \tau_{\text{shear}}$) driving fast-slow separation |
| $l_0$ | $1.0\,\mathrm{m}$ | Master turbulent mixing length scale |
| $N^2$ | $0.1\,\mathrm{s^{-2}}$ | Stratification: Brunt-Väisälä frequency squared |
| $\beta$ | $5.0$ | Stability coefficient governing buoyant suppression modification |
| $B_{0,\max}$ | $0.05\,\mathrm{m^2\,s^{-3}}$ | Maximum buoyant TKE destruction rate |
| $G_0$ | $0.3\,\mathrm{s^{-2}}$ | Synoptic/geostrophic wind shear forcing |
| $\gamma_s$ | $1.8\,\mathrm{m^{-2}\,s}$ | Turbulent momentum flux coefficient driving shear reduction |
| $\delta_{\text{reg}}, \eta_S$ | $10^{-2}, 10^{-4}$ | Numerical regularizations ensuring smooth $C^\infty$ manifolds |

---

**GSPT Manifold Geometry & Bifurcations**

* **Critical Manifold ($\mathcal{C}_0$):** Defined by $F(e,S) = 0$. Represents the quasi-equilibrium TKE surface that fast turbulence relaxes toward for any frozen shear value $S$.
* **Fast Stability & Hyperbolicity:** Governed by the fast eigenvalue $\lambda_{\text{fast}} = \frac{1}{\epsilon} \frac{\partial F}{\partial e}$.
* **Attracting Branch ($\mathcal{C}_0^+$):** $\lambda_{\text{fast}} < 0$. Perturbations in TKE decay rapidly; turbulence is dynamically stable.
* **Repelling Branch ($\mathcal{C}_0^-$):** $\lambda_{\text{fast}} > 0$. Perturbations grow; the equilibrium TKE state is unstable.


* **Saddle-Node Fold ($S_{\text{fold}}, e_{\text{fold}}$):** Identified by $F(e,S) = 0$ and $\frac{\partial F}{\partial e} = 0$. At this geometric boundary, normal hyperbolicity is lost ($\lambda_{\text{fast}} = 0$), and the restoring capability of turbulence vanishes.

---

**Key Meteorological Takeaways**

* **Manifold Bifurcation vs. $Ri_c$:** Regime transitions in the SBL are governed by loss of normal hyperbolicity along a continuous equilibrium manifold ($F = 0, F_e = 0$), rather than reaching a universal scalar Richardson number ($Ri_c \approx 0.25$).
* **Critical Slowing Down:** As an air mass approaches $S_{\text{fold}}$, $\lambda_{\text{fast}} \to 0^{-}$. High-frequency flux-tower data will exhibit increased variance and longer autocorrelation times in TKE prior to turbulence collapse.
* **Shear-TKE Feedback Loop:** Turbulent momentum mixing ($\gamma_s e S$) acts as a non-linear brake. High TKE reduces shear, pushing the trajectory leftward toward $S_{\text{fold}}$, triggering sudden transitions to weak or intermittent turbulence regimes.