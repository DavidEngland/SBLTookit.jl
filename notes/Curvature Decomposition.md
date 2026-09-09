### Curvature Decomposition

The central formulation of GSPT partitions the total observed physical vertical spatial curvature of a stability profile (such as the gradient Richardson number \(Ri_g(z)\)) into thermodynamic closure curvature and flux-coordinate mapping geometry.

* **Similarity Mapping:**
    \[\zeta(z) = \frac{z}{L(z)} \tag{88}\]
* **Coordinate Jacobian (\(\zeta_z \equiv \frac{\partial \zeta}{\partial z}\)):**
    \[\zeta_z = \frac{L(z) - z L'(z)}{L(z)^2} = \frac{1 - \zeta L'(z)}{L(z)} \tag{89, 110}\]
* **Coordinate Curvature (\(\zeta_{zz} \equiv \frac{\partial^2 \zeta}{\partial z^2}\)):**
    \[\zeta_{zz} = -\frac{z L''(z)}{L(z)^2} - \frac{2 L'(z)}{L(z)} \zeta_z = -\frac{L''(z)}{L(z)} \zeta - \frac{2 L'(z)}{L(z)} \zeta_z \tag{90, 110}\]
* **The Curvature Decomposition Equation:**
    \[\frac{d^2 Ri_g}{dz^2} = \underbrace{R''(\zeta) \zeta_z^2}_{C_{\text{constitutive}}} + \underbrace{R'(\zeta) \zeta_{zz}}_{C_{\text{mapping}}} + \mathcal{E}_{\Delta z} \tag{1, 88, 121}\]
  * _Where \(\mathcal{E}_{\Delta z}\) represents the explicit discrete accounting residual isolating grid discretization phase errors and non-commutation truncation._
* **The Linear Obukhov Profile Test (Transversality Gate):**
    If the vertical flux divergence is linear (\(L(z) = L_0 + az\)), no finite coordinate fold can exist:
    \[L''(z) = 0 \implies \zeta_z = \frac{L_0}{(L_0 + az)^2} > 0 \quad \forall z \in [0, \infty) \tag{120}\]
* **Nondegenerate Coordinate Fold Condition:**
    \[
        \zeta_z(z^_) = 0 \iff L(z^_) = z^{_}L'(z^_) \quad \text{under} \quad L''(z^_) \neq 0 \tag{8, 9, 90}
        \]
    At the fold locus \(z^_\), the total spatial curvature collapses to pure coordinate-stretching geometry:
    \[
        Ri_{zz}(z^_) = R'(\zeta^_) \zeta_{zz}(z^_) = -\frac{z^_ L''(z^_)}{L(z^_)^2 (1 + \beta_m \zeta^*)^2} \tag{11, 91}
\]
* **GSPT Fold Ratio Metric:**
    \[\text{Fold Ratio} = \frac{|C_{\text{mapping}}(z)|}{|C_{\text{constitutive}}(z)| + |C_{\text{mapping}}(z)| + \epsilon_c K_0} \tag{125, 131}\]

---

### Fast-Slow TKE–Shear Dynamical System

The SBL is modeled in 2D phase-space as a coupled singular perturbation system where microscale Turbulent Kinetic Energy (\(E\), the fast variable) relaxes rapidly toward the critical slow manifold sheet, while mesoscale vertical shear (\(S\), the slow variable) evolves slowly under geostrophic forcing.

* **Fast TKE Equation (\(E \ge 0\)):**
    \[\epsilon \frac{dE}{dt} = \mathcal{P} - \mathcal{B}_{bg}(E) - \mathcal{B}_{\text{mod}}(E) - \varepsilon_{\text{diss}} \tag{1, 599}\]
    \[\epsilon \frac{dE}{dt} = l_0 \sqrt{E + \delta} S^2 - l_0 \sqrt{E + \delta} \phi N^2 - l_0 \sqrt{E + \delta} c_b N^2 \frac{\alpha E}{(E + \alpha)^2} - \frac{E^{3/2}}{l_0} \tag{3, 5, 6, 8}\]
  * _Where \(\delta = 10^{-6}\text{ m}^2\text{s}^{-2}\) is the velocity-scale regularization floor restoring \(C^1\)-continuity of the fast vector field across the zero-energy limit._
* **Slow Wind Shear Equation (\(S \ge 0\)):**
    \[\frac{dS}{dt} = G - c_1 E S - c_2 S \tag{2, 600}\]
* **Fast Subsystem Stability Eigenvalue (\(\lambda_f\)):**
    \[\lambda_f(E, S) \equiv \frac{\partial f}{\partial E} = l_0 S^2 - \frac{2 B_{0,\max} E \delta_{\text{reg}}^2}{(E^2 + \delta_{\text{reg}}^2)^2} - \frac{3 \sqrt{E}}{2 l_0} \tag{95, 111}\]
* **Saddle-Node Fold Boundary of the Unregularized System:**
    \[x^3 + a x^2 - C = 0 \quad \left(\text{where } x \equiv S_{\text{fold}}^2, \ a \equiv \beta N^2, \ C \equiv \frac{27 B_{0,\max}^2}{4 l_0^4}\right) \tag{545, 557}\]
    \[e_{\text{fold}} = \frac{l_0}{\sqrt{3}} \sqrt{S_{\text{fold}}^2 + \beta N^2} = \frac{3 B_{0,\max}}{2 l_0 S_{\text{fold}}^2} \tag{545, 546}\]
* **Full 2D System Equilibrium Saddle-Node Condition (det \(J = 0\)):**
    \[F_e (\gamma e + r) - \gamma S F_S = 0 \tag{548}\]
* **Manifold Adhesion Distance Metric (\(d_\perp\)):**
    \[d_\perp(t) = \frac{|F(e,S)|}{\|\nabla_{e,S} F\| + \epsilon_{\text{floor}}} \tag{115, 148}\]

---

### SCM Implicit Coupling & Surface Energy Balance (SEB)

To avoid stiff numerical instabilities during rapid diurnal cooling, the prognostic Surface Energy Balance is analytically coupled directly into the tridiagonal matrix of the Backward-Euler SCM vertical diffusion stepper.

* **Continuous Surface Energy Balance:**
    \[C_s \frac{\partial T_s}{\partial t} = R_n(t) - H_0(t) - G_s(t) \tag{1, 177}\]
* **Boundary Heat and Soil Conduction Fluxes:**
    \[H_0 = \rho c_p K_h(z_1) \left( \frac{\theta_1 - T_s}{z_1} \right), \quad G_s = K_s (T_s - T_{\text{deep}}) \tag{177}\]
* **Analytical Elimination of \(T_s^{n+1}\) in Backward-Euler Step:**
    \[\gamma \equiv \frac{\rho c_p K_h^{n+1}}{z_1}, \quad \beta \equiv \frac{C_s}{\Delta t} + \gamma + K_s \tag{18}\]
    \[T_s^{n+1} = \frac{1}{\beta} \left[ \frac{C_s}{\Delta t} T_s^n + R_n^{n+1} + \gamma \theta_1^{n+1} + K_s T_{\text{deep}} \right] \tag{19}\]
  * _This \(T_s^{n+1}\) term is substituted directly into the atmospheric vertical diffusion row at level \(1\) (\(z_1\)), resolving surface coupling simultaneously without explicit stiffness._

---

### Monin–Obukhov Similarity and Profile Stability Functions

* **GLGS (eMOST / Modified Grachev et al. 2007) Stability Functions:**
    \[\phi_m(\zeta) = 1 + \frac{a_m \zeta}{(1 + b_m \zeta)^{2/3}} \tag{32, 419}\]
    \[\phi_h(\zeta) = Pr_0 \left(1 + \frac{a_h \zeta}{1 + b_h \zeta}\right) \tag{33, 419}\]
  * _With optimized GABLS3 campaign constants: \(Pr_0 = 0.98\), \(a_m = 5.0\), \(b_m = 0.3\), \(a_h = 5.0\), \(b_h = 0.4\)._
* **Integrated Stability Correction Functions (\(\psi_m, \psi_h\)):**
    \[\psi_m(\zeta) = -3 \frac{a_m}{b_m} \left[ (1 + b_m \zeta)^{1/3} - 1 \right] \tag{34, 420}\]
    \[\psi_h(\zeta) = -Pr_0 \frac{a_h}{b_h} \ln(1 + b_h \zeta) \tag{35, 420}\]
* **Symmetric SBL Quadratic Stability Functions (for \(Pr = 1\)):**
    \[f_m(Ri) = f_h(Ri) = \begin{cases} \left( \frac{Ri_c - Ri}{Ri_c} \right)^2, & 0 \le Ri \le Ri_c \ 0, & Ri > Ri_c \end{cases} \tag{31, 229}\]
  * _Where the quadratic derivative \(f_m'(Ri)\) remains continuous at \(Ri_c = 0.25\), preventing spatial and temporal numerical noise in first-order \(K\)-closure systems._

---

### Inverse-Problem Regularization & Gradient Extraction

To satisfy the **Tangential Cone Condition** and avoid noise amplification from non-linear quotients near weak-shear regions, GSPT applies **Track A Primitive Variable Regularization** directly to wind and temperature profiles prior to differentiation.

* **Wavenumber Transformation Coordinate:**
    \[\xi \equiv \ln(z / z_0) \tag{96, 562}\]
* **Tikhonov–Morozov Penalized Spline Functional:**
    \[\min_{f_s} \sum_{i=1}^{N_{\text{obs}}} \left( \frac{f_s(\xi_i) - q_i}{\sigma_i} \right)^2 + \alpha \int_0^{\xi_{\max}} [f_s''(\xi)]^2 d\xi \tag{96, 180}\]
* **Analytical Log-Height Derivatives:**
    \[\frac{df}{dz} = \frac{1}{z} \frac{df}{d\xi}, \quad \frac{d^2f}{dz^2} = \frac{1}{z^2} \left( \frac{d^2f}{d\xi^2} - \frac{df}{d\xi} \right) \tag{180, 562}\]
* **Tangential Cone Condition for Morozov Discrepancy Principle (MDP) Bisection Convergence:**
    \[\|F(x_2) - F(x_1) - F'(x_1)(x_2 - x_1)\|_Y \le \gamma \|F(x_2) - F(x_1)\|_Y \tag{159}\]

---

### Turbulence Quality-Control and Inversion Metrics

* **Non-Dimensional Triple-Point Dispersion (\(\delta_{\text{TP}}\)):**
    \[\delta_{\text{TP}}(t) = \frac{\max(z_K, z_e, z_{e_z}) - \min(z_K, z_e, z_{e_z})}{H_{\text{SBL}}(t)} \tag{331}\]
  * _Where \(z_K = \arg\min_z K_m(z, t)\), \(z_e(t) \implies e(z_e, t) = \alpha e_{\text{ref}}\), and \(z_{e_z} = \arg\max_z \left| \frac{\partial e}{\partial z} \right|\). Verified physical collapse events are gated where \(\delta_{\text{TP}} < 0.10\)._
* **Bulk Richardson (Sparse 2-Level Stencil Fallback):**
    \[Ri_b = \frac{g}{\overline{\theta_v}} \frac{(z_2 - z_1)(\theta_{v,2} - \theta_{v,1})}{(u_2 - u_1)^2 + (v_2 - v_1)^2} \tag{bulk logic}\]
* **Amplitude-Weighted CEOF Phase Gradient (Wave-Jet Separation):**
    \[\bar{\nabla \theta}_{\text{weighted}} = \frac{\sum_{i=1}^{N_z} |A(z_i)| \cdot \left| \frac{\partial \theta}{\partial z} \right|_{z_i}}{\sum_{i=1}^{N_z} |A(z_i)|} \tag{7, 196}\]
