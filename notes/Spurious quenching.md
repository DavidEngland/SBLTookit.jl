Spurious quenching occurs when non-linear stability functions interact with grid-scale gradient jumps, triggering an unphysical feedback loop that shuts down turbulent transport and causes runaway nocturnal cooling.

---

### 1. The Finite-Difference Amplification Mechanism

In numerical weather prediction (NWP) models, vertical mixing diffusivities $K_m$ and $K_h$ depend on local stability via empirical functions:

$$K_{m,h} = l^2 S_{m,h} \, f(Ri_g)$$

where $f(Ri_g)$ drops off steeply as $Ri_g$ increases toward a critical threshold ($Ri_{\text{crit}} \approx 0.25$).

Across adjacent grid levels $z_{k-1}$, $z_k$, and $z_{k+1}$ separated by grid spacing $\Delta z$, a Taylor expansion of $Ri_g$ reveals how spatial profile curvature drives discrete grid jumps:

$$Ri_g(z_k) = Ri_g(z_{k-1}) + \Delta z \left.\frac{\partial Ri_g}{\partial z}\right\vert{}_{k-1} + \frac{(\Delta z)^2}{2} Ri_{g,zz} + \mathcal{O}((\Delta z)^3)$$

At or near a coordinate fold ($\zeta_z \approx 0$), the first spatial derivative vanishes ($\frac{\partial Ri_g}{\partial z} \approx 0$), but the mapping curvature $Ri_{g,zz} = C_{\text{mapping}} = -R'(\zeta^*) \frac{z^* L''}{L^2}$ can be locally large.

Because $\Delta z$ in standard models is often 10–50 m in the lower SBL, the quadratic term $\frac{(\Delta z)^2}{2} Ri_{g,zz}$ causes $Ri_g$ to step abruptly over a single grid cell:

$$Ri_g(z_{k-1}) = 0.10 \quad (\text{turbulent}) \quad \longrightarrow \quad Ri_g(z_k) = 0.80 \quad (\text{quenched})$$

---

### 2. The Non-Linear "Clipping Cliff"

Most local planetary boundary layer (PBL) schemes (e.g., Mellor–Yamada, MYNN, Louis) employ stability functions $f(Ri_g)$ with extreme sensitivity near $Ri_{\text{crit}}$:

```
    f(Ri_g)
      1.0 |-------\
          |        \
          |         \  <-- Steep "Clipping Cliff"
      0.0 |----------\_________________ (K_min background floor)
          0        Ri_crit (~0.25)   Ri_g

```

When $Ri_g(z_k)$ is artificially pushed over $Ri_{\text{crit}}$ by $C_{\text{mapping}}$, $f(Ri_g)$ collapses from $\mathcal{O}(1)$ down to a tiny background diffusion floor $K_{\text{min}} \approx 10^{-4} \text{ m}^2/\text{s}$.

---

### 3. The Physical-Numerical Feedback Loop (Runaway Cooling)

Once $K_h(z_k)$ collapses to $K_{\text{min}}$, a self-amplifying numerical instability begins:

1. **Downward Heat Flux Shutdown:** Sensible heat flux $H_{k-1/2} = -\rho c_p K_h \left(\frac{\partial \theta}{\partial z}\right)_{k-1/2}$ drops to near zero across the interface between $z_{k-1}$ and $z_k$.
2. **Thermal Disconnection:** Ground longwave radiative cooling continues to remove heat from the lowest model layer ($z_{k-1}$). However, because interface $k-1/2$ is quenched, layer $k-1$ can no longer draw warm air down from layer $k$ above.
3. **Cascading Stabilization:** Layer $k-1$ cools rapidly, causing its local vertical potential temperature gradient $\frac{\partial \theta}{\partial z}$ to skyrocket.
4. **Downward Propagation:** The spike in $\frac{\partial \theta}{\partial z}$ pushes $Ri_g(z_{k-1})$ above $Ri_{\text{crit}}$ during the next time step, quenching layer $k-1$.

This creates a downward-propagating "extinction front" that locks the lowest grid cell into continuous, unphysical cooling, often overestimating surface cold biases by several degrees Kelvin.

---

### 4. Algorithmic Solutions

To prevent coordinate folds from triggering spurious quenching without removing physical stability limits, models can implement three targeted fixes:

* **Jacobian-Weighted Finite Differencing:** Evaluate discrete gradients over similarity space ($\Delta \zeta = \zeta_k - \zeta_{k-1}$) rather than physical space ($\Delta z$). Differencing across $\Delta \zeta$ absorbs $C_{\text{mapping}}$ directly, preventing artificial spikes in discrete $Ri_g$.
* **Curvature-Aware Flux Limiters:** Apply slope limiters (similar to Total Variation Diminishing schemes in fluid dynamics) to $Ri_g$ across cells where $\vert{}L(z) - z L'(z)\vert{} \to 0$, smoothing the profile before it enters $f(Ri_g)$.
* **$L_k$ Spatial Pre-Filtering:** Apply a 3-point vertical Shapiro filter to the discrete Obukhov length profile $L_k$ before computing stability functions, removing grid-scale noise in $L''_k$ while preserving the underlying SBL thermal structure.
