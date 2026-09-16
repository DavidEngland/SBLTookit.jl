### **Vertical Stratification & Layer Behavior**

* **Surface Layer (\(z \lesssim 0.1 z_i \sim 10\text{--}100\text{ m}\)):** Under classical Monin–Obukhov Similarity Theory (MOST), turbulent fluxes (\(u_*, H_0\)) are assumed nearly constant with height. In this shallow region, the surface Obukhov length \(L_0\) remains constant, making dimensionless stability \(\zeta(z) = z/L_0\) grow linearly with height. Under stable stratification (\(L_0 > 0\)), the gradient Richardson number \(Ri_g(z)\) increases monotonically until approaching the canonical critical threshold (\(Ri_c \approx 0.20\text{--}0.25\)).
* **Stable Outer Layer & Nocturnal Boundary Layer (\(z > 0.1 z_i\)):** Above the constant-flux surface layer, turbulent stress and sensible heat flux decay rapidly toward zero. As verified by field datasets (e.g., CASES-99, SHEBA, Cabauw/GABLS), surface MOST breaks down aloft. To describe this outer region, **Nieuwstadt’s local similarity theory** replaces surface parameters with a **local Obukhov length** \(L_{\text{loc}}(z) \equiv -\frac{u_*^3(z) \theta_v}{\kappa g \overline{w'\theta_v'}(z)}\) evaluated using local, \(z\)-dependent turbulent fluxes. Because local friction velocity \(u_*(z)\) decays faster with height than heat flux, \(L_{\text{loc}}(z)\) decreases aloft, causing local stability \(\zeta_{\text{loc}}(z) = z/L_{\text{loc}}(z)\) to grow nonlinearly.
* **Low-Level Jet (LLJ) Core (\(z \sim z_{\text{jet}}\)):** Near the nose of a nocturnal Low-Level Jet, vertical wind shear vanishes (\(S \equiv \left|\frac{\partial \mathbf{V}}{\partial z}\right| \to 0\)) while static temperature stratification remains positive (\(\theta_z > 0\)). This shear minimum causes the local gradient Richardson number to blow up toward infinity (\(Ri_g \to +\infty\)). In this region, surface-driven similarity scaling becomes physically meaningless, giving rise to an **"upside-down" boundary layer** where turbulent kinetic energy (TKE) is generated aloft by jet shear or gravity wave breaking and mixed downward toward the surface.
* **Convective Mixed Layer (\(L < 0\)):** During daytime heating, intense surface buoyancy drives convective plumes. Surface friction velocity \(u_*\) becomes secondary, and scaling transitions entirely from surface-layer Obukhov parameters to convective boundary layer scales (\(w_*\) and inversion height \(z_i\)).

---

### **Key Observational Lessons & NWP Model Implications**

1. **Local vs. Surface Similarity:** Field observations demonstrate that while surface similarity is restricted to the lowest \(10\%\) of the boundary layer, **local similarity theory (\(Ri = \mathcal{R}(\zeta_{\text{loc}})\)) accurately characterizes outer stable layers** provided that vertical profiles of turbulent stress \(u_*(z)\) and heat flux are predicted rather than assumed uniform.
2. **Supercritical Turbulence & Intermittency:** Early theoretical models assumed complete relaminarization and zero mixing whenever \(Ri_g > Ri_c \approx 0.25\). High-resolution field campaigns (such as SHEBA, VTMX, and TA-6) prove that **small-scale, non-Kolmogorov intermittent turbulence persists far above critical limits (\(Ri_g > 1.0\text{--}100\))**, driven by sub-grid shear, Kelvin–Helmholtz billows, internal gravity wave breaking, and radiative flux divergence.
3. **The Spurious Quenching Failure in NWP:** Enforcing hard surface-layer cutoffs (\(K_{m,h} = 0\) for \(Ri_g > Ri_c\)) across coarse numerical weather prediction (NWP) grids creates a severe numerical pathology. As under-resolved profile curvature artificially pushes discrete \(Ri_g\) across \(Ri_c\), unregularized closures drop eddy diffusivity to zero (\(K_h \to 0\)), shutting off downward heat delivery (\(H_0 \to 0\)). Deprived of heat from aloft, ground skin temperature \(T_s\) enters a downward spiral, driving unphysical surface cold-pool biases exceeding \(-3.5\text{ K}\) to \(-4.2\text{ K}\).

---

### **Mitigation via Upstream C-Z0HR Regularization**

To resolve this decoupling without applying ad-hoc empirical "long-tail" over-mixing curves, **Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR)** shifts the regularization operator upstream into the diagnostic scalar argument:

\[Ri_g^{\text{reg}} = Ri_c + (Ri_g - Ri_c) \left[ \frac{\Phi_\epsilon(Ri_g - Ri_c) + C(z, \Delta z)}{\Phi_\epsilon(Ri_g - Ri_c) + (1 + \alpha) C(z, \Delta z)} \right] \tag{1}\]

By absorbing grid-resolved spatial profile curvature \(C(z, \Delta z) = \frac{1}{2}\left|\frac{\partial^2 Ri_g}{\partial z^2}\right|(\Delta z)^2 + \epsilon_c Ri_c\) into \(Ri_g^{\text{reg}}\), C-Z0HR attenuates threshold sensitivity by \(\frac{1}{1+\alpha}\) (or \(66.7\%\) for \(\alpha = 2.0\)), capping discrete numerical loop gain (\(\mathcal{G} < 1\)) and preserving continuous surface-atmosphere thermal coupling while leaving baseline Monin–Obukhov similarity functions \(100\%\) untouched.

***

📉 **Would you like to run an SCM column diagnostic in `SBLToolkit.jl` comparing local Obukhov length profiles \(L_{\text{loc}}(z)\) against surface scaling \(L_0\) across a 24-hour GABLS3 diurnal cycle?**