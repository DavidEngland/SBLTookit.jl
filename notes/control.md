# Control Theory of Atmospheric Boundary Layers: Unified Lecture Notes

## Module 1: Fundamental Loop Gain & SCM Matrix Stability

The evolution of local stability across discrete time steps $n \to n+1$ in a Single-Column Model (SCM) is governed by a closed-loop feedback system. The continuous four-link physical gain chain evaluates to:

$$\mathcal{G}_{\text{physical}} = \underbrace{\left( \frac{g}{\theta_0 S^2} \right)}_{\text{Diagnostic Sensor}} \cdot \underbrace{\left( \frac{\partial \theta_z}{\partial H} \right)}_{\text{Plant Impedance}} \cdot \underbrace{\left( -\rho c_p \theta_z \right)}_{\text{Flux Actuator}} \cdot \underbrace{\left( l_h^2 S \, f_h'(Ri_g) \right)}_{\text{Closure Slope}}$$

This four-link scalar product represents a selected path gain. In a complete, discretized SCM, the heat flux perturbation depends on both diffusivity and potential temperature gradient:

$$\delta H = -\rho c_p \theta_z \, \delta K_h - \rho c_p K_h \, \delta \theta_z$$

* **Scalar Fixed-Point Stability ($\vert{}\mathcal{G}\vert{} < 1$):** Guarantees local contraction of perturbations at an isolated model interface.
* **Multidimensional Column Stability ($\rho(J) < 1$):** Dictates full SCM stability by requiring that the spectral radius of the state Jacobian $J = \frac{\partial \mathbf{x}^{n+1}}{\partial \mathbf{x}^n}$ remains strictly within the unit circle.

---

## Module 2: High-Gain Sensor Singularities & Actuator Switching

Numerical divergence in unregularized planetary boundary layer (PBL) schemes originates from the coupling of a singular sensor to a high-gain actuator switch:

* **LLJ Sensor Singularity:** At the nose of a Low-Level Jet where vertical wind shear vanishes ($S \equiv \left\vert{} \frac{\partial \mathbf{V}}{\partial z} \right\vert{} \to 0$), the diagnostic gain scales as $\mathcal{O}(S^{-2})$. Under constant shear-error variance, velocity profile noise is amplified as $\mathcal{O}(S^{-6})$, injecting high-frequency variance into the top of the gain chain.
* **Actuator Step-Discontinuities:** Monin–Obukhov Similarity Theory (MOST) closures derived by inverting similarity mappings with finite asymptotes $Ri_c \approx 0.20\text{--}0.25$ exhibit derivative divergence ($f_h'(Ri_g) \to -\infty$) as $Ri_g \to Ri_c^-$. The actuator acts as a hard off-switch, collapsing downward sensible heat flux ($H_0 \to 0$) and causing unphysical land-surface thermal decoupling and cold-pool runaway errors (biases exceeding $-3.5\text{ K}$ to $-4.2\text{ K}$).

---

## Module 3: Profile Curvature Kinematics & C-Z0HR Regularization

### Spatial Curvature Kinematics

Expanding $Ri(z) = \mathcal{R}(\zeta(z))$ with Monin–Obukhov variable $\zeta(z) = z / L(z)$ yields:

$$Ri_{zz} \equiv \frac{\partial^2 Ri}{\partial z^2} = \mathcal{R}''(\zeta) (\zeta_z)^2 + \mathcal{R}'(\zeta) \zeta_{zz}$$

In the constant-flux surface layer ($\zeta_z = 1/L$, $\zeta_{zz} = 0$), this simplifies to $Ri_{zz} = \frac{1}{L^2} \mathcal{R}''(\zeta)$. For standard MOST closures with $\mathcal{R}(\zeta) = \frac{\zeta}{1 + \beta \zeta}$, the second derivative $\mathcal{R}''(\zeta) = -\frac{2\beta}{(1 + \beta \zeta)^3}$ remains strictly negative and smooth.

### Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR)

To control feedback gain without altering empirical constants ($\beta = 5.0$), C-Z0HR inserts a partial slope attenuator $\mathcal{D}_\alpha^{(R)}(Ri_g) \equiv \left. \frac{\partial \widetilde{Ri}_g}{\partial Ri_g} \right\vert{}_C$ into the closure link, yielding $\mathcal{G}_{\text{reg}} = \mathcal{G}_{\text{physical}} \cdot \mathcal{D}_\alpha^{(R)}(Ri_g)$.

With grid curvature metric $C = \frac{1}{2}\vert{}Ri_{g,zz}\vert{}(\Delta z)^2 + \epsilon_c Ri_c$, the exact-threshold mapping slope is:

$$\mathcal{D}_\alpha^{(R)}(Ri_c) = \frac{\epsilon_c + C}{\epsilon_c + (1 + \alpha)C}$$

* **Coarse Grids / High Curvature ($C \gg \epsilon_c$):** $\mathcal{D}_\alpha^{(R)}(Ri_c) \approx \frac{1}{1 + \alpha}$. Setting $\alpha = 2.0$ attenuates threshold sensitivity by $66.7\%$ ($\mathcal{D}_2^{(R)} = 1/3$), enforcing contraction mapping stability.
* **Fine Grids ($\Delta z \to 0$):** Physical profile curvature vanishes ($C \to \epsilon_c Ri_c$), recovering exact physical similarity laws away from exact threshold.

---

## Module 4: Geometric Surface Pattern Theory (GSPT) & Fast–Slow Systems

Atmospheric boundary layer columns can be cast into a singularly perturbed fast–slow state space ($\epsilon \ll 1$):

$$\epsilon \frac{de}{dt} = \mathcal{P}(e, S) - \mathcal{B}(e) - \varepsilon_{\text{diss}}(e) \quad \text{(Fast Variable: TKE } e\text{)}$$

$$\frac{dS}{dt} = G_0 - \gamma_s e S - r_s S \quad \text{(Slow Variable: Shear } S\text{)}$$

* **Normal Hyperbolicity:** The critical manifold $\mathcal{C}_0$ (where $\mathcal{P} - \mathcal{B} - \varepsilon_{\text{diss}} = 0$) is attracting when the fast eigenvalue is negative:

$$\lambda_f(e, S) = \frac{1}{\epsilon} \frac{\partial F}{\partial e} < 0$$

* **The Decoupling Principle:** Spatial coordinate turning points ($\zeta_z = 0$) resulting from jet-nose coordinate compression are logically independent of fast eigenvalue zero-crossings ($\zeta_z = 0 \quad \not\!\!\!\!\!\iff \lambda_f = 0$). Gating relaminarization triggers on $\lambda_f \ge 0$ prevents SCMs from treating spatial coordinate geometry as physical turbulence collapse.

---

## Module 5: Multilayer Conditioning & 4D-Var Assimilation

| Regime | Mathematical Condition | Numerical / Physical Impact |
| --- | --- | --- |
| **Implicit SCM Solvers** | $\kappa(J) < \infty$, where $J = I - \Delta t \frac{\partial \mathbf{F}}{\partial \mathbf{x}}$ | Eliminates Newton–Raphson iteration stalls by bounding off-diagonal coupling terms ($\frac{\partial H_{k\pm 1/2}}{\partial \theta_{z,k}}$). |
| **4D-Var Data Assimilation** | Continuous $C^\infty$ operator smoothness | Eliminates gradient singularities at $Ri_c$, enabling stable Tangent Linear Model (TLM) and Adjoint propagation ($\frac{\partial J_{\text{cost}}}{\partial \mathbf{x}}$). |
