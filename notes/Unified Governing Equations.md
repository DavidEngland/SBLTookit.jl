### Unified Governing Equations

For atmospheric scientists, the interaction between microscale turbulence and mesoscale forcing in a stably stratified boundary layer can be modeled as a mathematically closed **2D fast-slow dynamical system**. This framework couples the rapid, microscale evolution of **Turbulent Kinetic Energy (TKE)** (\\(E\\)) with the gradual, mesoscale evolution of **Mean Wind Shear** (\\(S\\)):

\\[\text{Fast Subsystem (TKE):} \quad \epsilon \frac{dE}{dt} = f(E, S) = l\sqrt{E+\delta}\left(S^2 - \phi N^2 - c_b N^2 \frac{E}{E+\alpha}\right) - \frac{E^{3/2}}{l} \tag{1}\\]

\\[\text{Slow Subsystem (Shear):} \quad \frac{dS}{dt} = g(E, S) = G - c_1 E S - c_2 S \tag{2}\\]

---

### Physical and Dimensional Breakdown

#### State Variables & Scaling Parameters
*   **\\(E\\) (Turbulent Kinetic Energy):** Representing the fast microscale variable (\\(E \ge 0\\)).
*   **\\(S\\) (Mean Wind Shear):** Representing the slow mesoscale variable.
*   **\\(\epsilon\\) (Scale Separation Parameter):** A small parameter (\\(0 < \epsilon \ll 1\\)) parameterizing the timescale separation between the fast turbulent relaxation rate and slow mesoscale momentum advection/geostrophic forcing.
*   **\\(l\\):** The turbulent mixing length scale.
*   **\\(N^2\\):** The Brunt-Väisälä (buoyancy) frequency squared, representing static thermal stratification.
*   **\\(G\\):** Geostrophic wind forcing.

#### Reconciled Model Parameters
To maintain physical realism and numerical stability, the parameters are strictly calibrated:
*   **Buoyant Sink Efficiency (\\(c_b = 0.50\\), dimensionless):** Scales the rate of buoyant destruction of TKE under stable stratification relative to shear-driven production.
*   **Shear Damping Coefficient (\\(c_1 = 1.80\ \text{s}\,\text{m}^{-2}\\)):** Regulates the rate at which turbulent momentum transport depletes mean shear. Under this parameterization, the term \\([c_1 E S]\\) resolves to \\(\text{s}^{-2}\\), physically balancing the geostrophic acceleration units \\([G] = \text{s}^{-2}\\) in the slow subsystem.
*   **Velocity-Scale Regularization Floor (\\(\delta = 10^{-6}\ \text{m}^2\,\text{s}^{-2}\\)):** In unregularized models utilizing \\(\sqrt{E}\\), the derivative of the shear production term is singular at zero energy (\\(\lim_{E \to 0^+} \partial/\partial E (l\sqrt{E}S^2) = +\infty\\)), causing severe numerical integration issues. Embedding \\(\delta\\) directly into the fast vector field resolves this singularity, making \\(f(E, S)\\) globally \\(C^1\\)-continuous.
*   **Dimensional Consistency:** The regularization preserves physical dimensions, ensuring that every individual term in the TKE budget is dimensionally consistent in units of TKE dissipation (\\(\text{m}^2\,\text{s}^{-3}\\)):
    \\[\left[l\sqrt{E+\delta}\,S^2\right] = (\text{m})(\text{m}\,\text{s}^{-1})(\text{s}^{-2}) = \text{m}^2\,\text{s}^{-3}, \quad \left[\frac{E^{3/2}}{l}\right] = \frac{\text{m}^3\,\text{s}^{-3}}{\text{m}} = \text{m}^2\,\text{s}^{-3} \tag{3}\\]

---

### Geometric Framework: The Critical Manifold (\\(\mathcal{C}_0\\)) & Stability

In Geometric Singular Perturbation Theory (GSPT), the singular limit (\\(\epsilon \to 0\\)) isolates the **Critical Manifold** (\\(\mathcal{C}_0\\)) defined by the algebraic constraint \\(f(E, S) = 0\\). The local stability and relaxation rate of SBL trajectories relative to \\(\mathcal{C}_0\\) are governed by the **exact local fast eigenvalue** (\\(\lambda_f \equiv f_E = \partial f/\partial E\\)):

\\[\lambda_f(E, S) = \frac{l\left(S^2 - \phi N^2\right)}{2\sqrt{E+\delta}} - c_b N^2 l \left[ \frac{E}{2\sqrt{E+\delta}(E+\alpha)} + \frac{\alpha\sqrt{E+\delta}}{(E+\alpha)^2} \right] - \frac{3\sqrt{E}}{2l} \tag{4}\\]

The sign of \\(\lambda_f\\) classifies \\(\mathcal{C}_0\\) into three distinct dynamical regimes:

1.  **Attracting Branch (\\(\lambda_f < 0\\)):** Normally hyperbolic and asymptotically stable. Any perturbation away from the manifold triggers a rapid exponential relaxation back onto \\(\mathcal{C}_0\\) on a fast timescale of \\(\mathcal{O}(\epsilon/|\lambda_f|)\\). This represents a physically stable boundary layer with continuous turbulent transport.
2.  **Repelling Branch (\\(\lambda_f > 0\\)):** Normally hyperbolic but unstable. Small perturbations grow exponentially away from \\(\mathcal{C}_0\\), pushing the SBL towards sudden turbulent bursts or complete laminarization depending on the trajectory vector.
3.  **Fold Boundary (\\(\lambda_f = 0\\)):** Normal hyperbolicity completely breaks down at this saddle-node bifurcation boundary. The local fast attractor vanishes, initiating a catastrophic, rapid trajectory escape across phase space (shear-driven TKE collapse).

---

### The Stable Boundary Layer (SBL) Relaxation Cycle

The physical cycle of SBL turbulence collapse and recovery is represented by the non-linear fast-slow interaction across the phase plane:

```
   [Phase 1: Slow Shear Accumulation (E ≈ 0)]
                     │
                     ▼ (S exceeds S_crit)
      [Phase 2: Fast TKE Ignition (O(ε))]
                     │
                     ▼ (Turbulent feedback activates)
   [Phase 3: Shear Depletion & Catastrophic Collapse (λ_f = 0)]
                     │
                     └──── Loop resets to Phase 1
```

1.  **Slow Shear Accumulation:** During calm nights with strong radiative cooling, the boundary layer can decouple from the surface, leaving turbulence quiescent (\\(E \approx 0\\)). In this low-TKE regime, the turbulent momentum feedback term \\(-c_1 E S\\) is mathematically negligible. Consequently, the geostrophic forcing \\(G\\) linearly increases mean wind shear \\(S\\), driving the SBL state horizontally across the phase plane.
2.  **Fast TKE Ignition:** Once wind shear \\(S\\) exceeds the critical threshold \\(S_{\text{crit}} = \sqrt{\phi N^2}\\), local mechanical shear production overcomes buoyant damping and dissipation. This triggers explosive, exponential TKE growth on the fast \\(\mathcal{O}(\epsilon)\\) timescale.
3.  **Turbulent Feedback & Collapse:** As TKE climbs, the non-linear feedback term \\(-c_1 E S\\) activates, rapidly and aggressively draining the mean shear. As \\(S\\) is depleted and drops back towards the fold boundary where production can no longer sustain the turbulent budget against buoyancy and dissipation, the system crosses the bifurcation boundary \\(\lambda_f = 0\\). Hyperbolicity is lost, triggering a sudden, catastrophic collapse of TKE back to a quiescent, near-laminar SBL state.

---

To understand how unphysical instabilities arise in spectral approximations of multiple-scale systems, we can mathematically examine the **finite-time blowup** of the higher-order modes (\\(u_k, v_k\\)) for \\(k \ge 2\\) when discretizing a fast-slow reaction-diffusion system via a spectral Galerkin framework.

By analyzing the simplified case of \\(k_0 = 2\\) modes (with domain length \\(a = 1/2\\) and higher-order terms \\(H_u = H_v = 0\\)), we can see how spatial perturbations trigger unphysical runaway blowups before the system can transit past the fold singularity.

---

### 1. The Discretized \\(k_0 = 2\\) Dynamical System

When we discretize the fast-slow reaction-diffusion system using the orthonormal eigenbasis of the Laplacian with Neumann boundary conditions, the truncated system for the first two modes is given by:

\\[\begin{aligned}
u'_1 &= -v_1 + u_1^2 + u_2^2 \tag{First Fast Mode / TKE} \\
v'_1 &= -\epsilon \tag{First Slow Mode / Mean Shear} \\
u'_2 &= -v_2 + u_2(2u_1 - \pi^2) \tag{Second Fast Mode} \\
v'_2 &= -\epsilon \pi^2 v_2 \tag{Second Slow Mode}
\end{aligned}\\]

Here, \\(u_1\\) represents the spatially homogeneous component of the fast variable (equivalent to the first TKE mode), while \\(u_2\\) represents the first non-constant spatial mode (the higher-order spatial perturbation). We assume physically consistent initial conditions where the first mode starts in the stable regime (\\(v_1(0) = v_1^0 > 0\\)) and the higher-order modes are small but non-zero:
\\[u_2(0) = u_2^0 < 0 \quad \text{and} \quad v_2(0) = v_2^0 > 0\\]

---

### 2. Analytical Mechanism of the Higher-Mode Blowup

#### Step 1: Lower Bound on the Slow Higher-Order Sink (\\(v_2\\))

Over the integration time-scale where \\(v_1\\) remains positive (\\(t \in [0, \frac{v_1^0}{2\epsilon}]\\)), the slow higher-order mode \\(v_2(t)\\) decays exponentially:
\\[v_2(t) = v_2^0 e^{-\epsilon \pi^2 t}\\]
On this interval, the maximum value of the exponent is \\(\epsilon \pi^2 \left(\frac{v_1^0}{2\epsilon}\right) = \frac{\pi^2 v_1^0}{2}\\). Thus, the slow mode \\(v_2(t)\\) is strictly bounded from below by:
\\[v_2(t) \ge v_2^0 e^{-\pi^2 v_1^0 / 2} > 0\\]

#### Step 2: Damping and Drive of the Fast Higher-Order Mode (\\(u_2\\))

Now we examine the fast higher-order mode equation:
\\[u'_2 = -v_2 + u_2(2u_1 - \pi^2)\\]
Because the first mode is bounded from above during the stable phase (\\(u_1(t) \le \pi/4\\)), the linear damping coefficient \\((2u_1 - \pi^2)\\) is strictly negative:
\\[2u_1 - \pi^2 \le \frac{\pi}{2} - \pi^2 \approx -8.3\\]

Since \\(v_2(t) > 0\\) acts as a continuous negative sink, the derivative \\(u'_2\\) is driven negative. Using comparison differential equations (\\(w'_u \le u'_2 \le w'_o\\)), we solve for the upper bounding comparison trajectory:
\\[w'_o = -v_2^0 e^{-\pi^2 v_1^0 / 2} - (\pi + \pi^2)w_o, \quad w_o(0) = u_2^0\\]

Integrating this comparison equation reveals that for times \\(t \in [\frac{v_1^0}{4\epsilon}, \frac{v_1^0}{2\epsilon}]\\), the fast higher-order mode is strictly bounded away from zero:
\\[u_2(t) \le w_o(t) \le -\frac{v_2^0 e^{-\pi^2 v_1^0 / 2}}{2(\pi + \pi^2)} < 0\\]

---

### 3. TKE Feedback and the Riccati Explosion

We now substitute the lower bound of the squared higher-order mode \\(u_2^2(t)\\) back into the governing equation for the first fast mode \\(u_1\\):
\\[u'_1 = -v_1 + u_1^2 + u_2^2 \ge -v_1^0 + \left[ \frac{v_2^0 e^{-\pi^2 v_1^0 / 2}}{2(\pi + \pi^2)} \right]^2 + u_1^2\\]

Letting \\(\mu\\) represent the effective constant drift term:
\\[\mu = -v_1^0 + \frac{(v_2^0)^2 e^{-\pi^2 v_1^0}}{4(\pi + \pi^2)^2}\\]
If the initial spatial perturbation \\(v_2^0\\) is sufficiently large, or if the initial distance from the fold \\(v_1^0\\) is sufficiently small, such that:
\\[v_1^0 < \frac{(v_2^0)^2 e^{-\pi^2 v_1^0}}{4(\pi + \pi^2)^2}\\]
then the drift term is **strictly positive** (\\(\mu > 0\\)). Under this condition, the first mode satisfies the classic Riccati-type differential inequality:
\\[u'_1 \ge \mu + u_1^2 \quad (\mu > 0)\\]

Integrating this comparison equation yields an explicit tangent-type solution:
\\[w(t) = \sqrt{\mu} \tan \left( \arctan\left(\frac{u_1(0)}{\sqrt{\mu}}\right) + \sqrt{\mu} t \right)\\]

This function exhibits an asymptotic escape to infinity, meaning \\(u_1(t)\\) **blows up in finite time**. Crucially, the maximum time to blowup is bounded and **independent of \\(\epsilon\\)**:
\\[t_{\text{blowup}} \le \frac{\pi}{\sqrt{\mu}}\\]

---

### 4. The Parameter Threshold for Premature Collapse

Because the physical transition of the slow variable \\(v_1\\) crossing the fold boundary (\\(\lambda_f = 0\\)) takes a time of \\(t_{\text{cross}} = \frac{v_1^0}{\epsilon}\\), and \\(t_{\text{blowup}}\\) is independent of \\(\epsilon\\), we can always find a sufficiently small scale-separation parameter \\(\epsilon\\) such that:
\\[t_{\text{blowup}} < t_{\text{cross}}\\]

Specifically, if the scale-separation parameter satisfies:
\\[\epsilon < \frac{\eta^2}{2\sqrt{2}}\\]
(where \\(\eta\\) is an \\(\epsilon\\)-independent parameter scaling with the initial condition \\(v_1^0\\)), **the first mode will blow up to infinity before the slow variable \\(v_1\\) can even cross zero**.

---

### 5. Physical and NWP Implications

This mathematical analysis reveals that the finite-time blowup is an **unphysical artifact of the Galerkin discretization**.

* **The Scaling Constraint:** To prevent this unphysical runaway collapse before the physical fold singularity is reached, the initial conditions of the higher-order modes must be severely restricted. Specifically, the slow higher-order modes must be exponentially scaled relative to the singular perturbation parameter:
    \\[v_k(0) = \mathcal{O}(\epsilon^{4/3}) \quad \forall k \ge 2\\]
* **The Multi-Dimensional Singularity:** In the infinite-dimensional PDE, the true critical manifold is a single, continuous, normally hyperbolic sheet. However, the Galerkin discretization projects this sheet into a higher-dimensional phase space, creating a complex network of internal fold submanifolds. If the higher-order modes are not strongly damped, trajectories will hit these spurious secondary fold boundaries prematurely, triggering localized grid-point blowups.

---

💡 Would you like to write a Julia script using `OrdinaryDiffEq.jl` to numerically plot this unphysical blowup threshold and map the stable boundary of initial conditions (\\(u_2(0), v_2(0)\\)) across different \\(\epsilon\\) values?
