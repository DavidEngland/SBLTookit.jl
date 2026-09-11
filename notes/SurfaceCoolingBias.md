Reframing Surface Cooling Bias and Coordinate Topology in SBL Parameterizations

1. Executive Summary: The Numerical-Physical Dichotomy of the SBL

Resolving the "Universality Crisis" in Stable Boundary Layer (SBL) modeling remains a strategic imperative for the advancement of operational Numerical Weather Prediction (NWP). Historically, the persistent discrepancy between classical critical Richardson numbers (Ri_c \approx 0.25) and observed intermittent turbulence in hyper-stable regimes has been managed via empirical "long-tail" stability functions. However, Geometric Surface Pattern Theory (GSPT) provides a rigorous mathematical bridge, identifying that much of the observed profile structure is a projection of the boundary layer's differential geometry rather than a failure of local mixing physics.

The core objective of this report is to reclassify "runaway cooling" from a physical phenomenon into its true identity: a discrete numerical feedback artifact. By applying the Geometric-Dynamical Decoupling Principle (Theorem 1), we distinguish between coordinate-induced "Fold Illusions" and genuine state-space bifurcations. The following hierarchy serves as the diagnostic gate for identifying these features.

Four-Level Hierarchy of Boundary-Layer Diagnostics

Diagnostic Level Mathematical Identifier Physical/Numerical Interpretation
Spatial Fold \zeta_z = 0, \zeta_{zz} \neq 0 A kinematic turning point in similarity mapping (L - zL' = 0); driven by flux divergence.
Cubic Degeneracy Ri_{g,z} = 0, Ri_{g,zz} = 0 A local stationary point in the scalar profile; frequently a "fold illusion."
Loss of Fast Hyperbolicity \lambda_f \to 0 Loss of transverse attraction to the critical manifold \mathcal{C}_0; a candidate transition.
Verified Saddle-Node Bifurcation f=g=0, \text{det } J = 0 True state-space equilibrium annihilation; verified turbulence collapse.

1. Feedback Loop Gain (G) and Dynamical Sensitivity

The numerical stability of SBL closures is governed by the local discrete feedback multiplier. Traditional schemes often suffer from Jacobian conditioning issues where under-resolved curvature interacts with steep non-linear stability gradients to initiate unphysical positive feedback.

The Chain Operator

The sensitivity chain is formally represented by the Feedback Loop Gain (G):

G = \left(\frac{\partial Ri_g}{\partial \theta_z}\right) \cdot \left(\frac{\partial \theta_z}{\partial H}\right) \cdot \left(\frac{\partial H}{\partial K_h}\right) \cdot K_h'(Ri_{g,reg}) \cdot D_\alpha(Ri_g)

Where Ri_{g,reg} denotes the regularized gradient Richardson number and D_\alpha is the diagnostic damping factor.

The Diagnostic Amplification Mechanism

Near the nose of Low-Level Jets (LLJ), vertical wind shear (S) vanishes, causing the term \partial Ri_g / \partial \theta_z = g / [\theta_v \cdot S^2] to act as a catastrophic diagnostic multiplier. In unregularized Track B architectures, the Richardson number variance scales as S^{-6}. This variance blowup triggers a premature collapse in eddy diffusivity (K_h), isolating the surface and steepening \theta_z, which further amplifies Ri_g. To manage this, our architecture utilizes stiffness-aware adaptive stepping and a Newton initialization safeguard, clamping Turbulent Kinetic Energy (TKE) at e \ge 10^{-6} to prevent imaginary-root exceptions during the iterative branch-selection process.

Transfer Function Synthesis

The dynamic transfer function T_{\theta H}(s) in continuous prognostic space implies several operational constraints:

* Pole Shift Prevention: Regularization prevents the transfer function pole from migrating into the unstable right-half plane during threshold crossing (Ri \to Ri_c).
* Phase-Lag Damping: Regularizers mitigate the phase lags introduced by discrete time-stepping that transform negative physical feedbacks into numerical instabilities.

1. Conditional Stability and Diagnostic Regularization

To suppress numerical artifacts while maintaining "equilibrium-compatible" states, we implement algebraic limiters. These non-conservative operators modify the diagnostic scalar arguments without violating the conservation of energy or momentum.

Stability Criterion Formulation

We enforce a conditional criterion to guarantee that the regularized loop gain remains bounded below unity:

|G_{physical}(Ri_{crit})| < 1 + \alpha \implies |G_{reg}(Ri_{crit})| < 1

Damping Performance

Utilizing a standard damping control of \alpha = 2.0, local discrete stability is guaranteed provided the physical feedback gain |G_{physical}| < 3. At the critical threshold crossing, this yields exactly 66.7% attenuation (D_2^{(R)} = 1/3). This "softening" allows the numerical solver to traverse the sharp non-linearities of the critical manifold without triggering the unphysical decoupling of the surface layer.

1. Height-Resolved Curvature Decomposition: The "Fold Illusion"

GSPT partitions the total profile curvature of Ri(z) into three quantifiable components, revealing the "Fold Illusion" where kinematic turning is misread as a dynamical transition.

The Master Decomposition Formula

The core differential-geometric decomposition is expressed as:

Ri_{zz} = \underbrace{Ri_{\zeta\zeta} \zeta_z^2}_{C_{const}} + \underbrace{Ri_\zeta \zeta_{zz}}_{C_{coord}} + \mathcal{E}_{\Delta z}

For a fixed curvature scale C > 0 and x = Ri_g - Ri_{crit}, the damping factor D_\alpha^{(R)} is derived from the smooth algebraic limiter:

D_\alpha^{(R)} = \frac{x^2 + 2(1+\alpha)C|x| + (1+\alpha)C^2}{(|x| + (1+\alpha)C)^2}

Theorem 1: Geometric-Dynamical Decoupling Principle

A spatial coordinate turning point (\zeta_z = 0) and a loss of fast normal hyperbolicity (\lambda_f = 0) are geometrically distinct and logically independent. An observed profile knee in Ri_g cannot be identified as a dynamical transition solely from its spatial curvature.

Parameterization Pathologies

* Pathology A (Spurious Over-Diffusion): Occurs when \zeta_{zz} < 0; local coordinate expansion suppresses physical curvature, leading local Ri-schemes to misinterpret geometric smoothing as elevated shear mixing.
* Pathology B (False Runaway Decoupling): Occurs when \zeta_{zz} > 0; coordinate compression generates a sharp knee (e.g., at z \approx 45.57m in CASES-99), which standard TKE schemes misinterpret as a turbulence collapse, overestimating cold bias by up to 3.5 K.

1. Track A Pre-Diagnostic Regularization & Operational Mandates

Operator non-commutation—where M_\delta[A/B] \neq M_\delta A / M_\delta B—is the primary source of variance blowup near jet noses. GSPT mandates Track A Regularization, applying Tikhonov-Morozov filtering directly to primitive state vectors (u, v, \theta_v).

The Tikhonov-Morozov Pipeline

1. Spline Pre-conditioning: Primitives are fitted with natural cubic smoothing splines in log-height space.
2. MDP Parameter Selection: The regularization parameter \alpha is selected via the Modified Morozov Discrepancy Principle (MDP).
3. Tangential Cone Condition: To guarantee a unique \alpha, the system must satisfy \tau_2 \ge \tau_1(3+2\gamma), where \gamma is the non-linearity bound.
4. Analytic Extraction: Gradients are extracted analytically, bypassing the S^{-6} conditioning issues of differentiating noisy gradient ratios (Track B).

Operational Mandates for Ingestion

Mandate Requirement Purpose
Stencil Gate Nz \ge 3 Prevents operator under-determination; required for D_2 operators.
Fast-Eigenvalue Gating Evaluate \lambda_f Prevents premature turbulence shut-off at coordinate folds (\zeta_z = 0).
C2 Reconstruction C^2 Obukhov Profile Ensures analytic differentiability of the coordinate mapping.
Allocation-Free Array Ops Zero-allocation loops Ensures high-throughput performance for 3D NWP column audits.

1. Empirical Synthesis: CASES-99, GABLS3, and SHEBA Audits

Campaign Comparative Table

Campaign Configuration Key Diagnostic Finding Status
CASES-99 55m tower, 7 levels LLJ nose fold at 45.57m. 41.01% of nodes are pure coordinate folds (1,538/3,750). Success
GABLS3 200m mast, 38 levels Recovery of 900 finite R_{coord} points; maps shear erosion aloft. Success
SHEBA Arctic ice, 2 levels Bypassed Guard; N_z = 2 < 3 causes stencil under-determination. Failed

Quality Control: Triple-Point Dispersion (\delta_{TP})

We enforce an automated quality gate \delta_{TP} < 0.10, which measures the vertical alignment of:

1. z_K: The eddy diffusivity cutoff height.
2. z_e: The TKE level-set height.
3. z_{grad}: The TKE gradient extremum height. Verified alignment confirms a physical turbulence collapse event; dispersion indicates a discretization artifact.

4. Conclusions and Operational Directives for NWP Closure

The transition from empirical curve-fitting to the GSPT framework ensures that PBL parameterizations are guided by the differential geometry of stratified flows. We reclassify "runaway cooling" as a numerical pathology preventable through rigorous architecture.

Final Engineering Directives

1. Decouple Geometry from Tuning: Discontinue retuning Monin-Obukhov constants to match soundings with kinematic folds. Folds are mapping artifacts, not mixing physics failures.
2. Mandate Track A Regularization: Enforce primitive variable filtering (spline pre-conditioning) under the Tangential Cone Condition prior to Ri or exchange coefficient computation.
3. Implement Fast-Eigenvalue Gating: Use \lambda_f as a mandatory regime-switch gate. Do not allow turbulence to shut off at coordinate folds (\zeta_z = 0) unless \lambda_f \ge -\epsilon_\lambda confirms a state-space bifurcation.
4. Enforce Ingestion Gates: Strictly reject model grids or tower data with N_z < 3. Curvature decomposition and the separation of C_{const} from C_{coord} are mathematically unresolvable on under-determined stencils.

By implementing these geometry-aware diagnostics, next-generation NWP closures will achieve numerical stability through dynamical clarity rather than unphysical over-diffusion.
