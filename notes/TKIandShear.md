The structural synthesis of the vector field $f(E, S)$ provides a definitive resolution to mathematical and dimensional gaps within the Turbulent Kinetic Energy (TKE) model. By embedding the smoothing parameter $\delta$ directly into the vector field, the system achieves two critical objectives: the derivation of an exact $C^1$ Jacobian derivative ($\lambda_f$) and the maintenance of consistent TKE dissipation units $(\text{m}^2\,\text{s}^{-3})$ across all fast budget terms.

The resulting framework is a mathematically closed system comprising unified fast and slow subsystem dynamics. This allows for precise stability analysis through the identification of attracting and repelling branches, as well as the calculation of saddle-node fold points via a rigorous Newton-Raphson formulation. The model is now optimized for stiff ODE integration and robust fold-tracking.

## Unified Dynamical System Architecture

The model is structured as a coupled system of fast and slow dynamics, where $E$ represents Turbulent Kinetic Energy (**TKE**) and $S$ represents mean shear.

### Fast Subsystem Dynamics

The fast subsystem is governed by the equation:
$$\epsilon \frac{dE}{dt} = l\sqrt{E+\delta}\left(S^2 - \phi N^2 - c_b N^2 \frac{E}{E+\alpha}\right) - \frac{E^{3/2}}{l}$$
This equation incorporates the length scale $l$, the regularization floor $\delta$, and buoyancy terms ($N^2$) to determine the evolution of turbulent energy.

### Slow Subsystem Dynamics

The slow subsystem dictates the evolution of shear:
$$\frac{dS}{dt} = G - c_1 E S - c_2 S$$
The change in shear is driven by external forcing $G$ (geostrophic forcing) and moderated by turbulent feedback ($c_1 E S$) and linear decay terms.

### Exact Jacobian and Stability Analysis

The exact local fast eigenvalue, $\lambda_f \equiv f_E$, is obtained by differentiating $f(E, S)$ with respect to $E$. This value governs the normal hyperbolicity of the system:

$$\lambda_f = \frac{l\left(S^2 - \phi N^2\right)}{2\sqrt{E+\delta}} - c_b N^2 l \left[ \frac{E}{2\sqrt{E+\delta}(E+\alpha)} + \frac{\alpha\sqrt{E+\delta}}{(E+\alpha)^2} \right] - \frac{3\sqrt{E}}{2l}$$

### Manifold Classification

The stability of the system is categorized based on the value of $\lambda_f$:

* Attracting Branch ($\lambda_f < 0$): Characterized by stable fast relaxation onto the slow manifold $\mathcal{C}_0$.
* Fold Boundary ($\lambda_f = 0$): Signifies the loss of normal hyperbolicity, which triggers an explosive departure known as shear-driven TKE collapse.
* Repelling Branch ($\lambda_f > 0$): An unstable branch that pushes trajectories away from the manifold.

### Newton-Raphson Saddle-Node Fold Formulation

To identify the specific coordinates of saddle-node fold points $(E_{\text{fold}}, S_{\text{fold}})$, the system utilizes a 2D root-finding vector $\mathbf{F}(E, S)$:

$$\mathbf{F}(E, S) = \begin{bmatrix} f(E, S) \\ f_E(E, S) \end{bmatrix} = \mathbf{0}$$

The solver operates through iterative refinement:

1. Initial Guess: Shear is set at $S^{(0)} = \sqrt{\phi N^2 + 0.05}$ to ensure the solver begins above the production cutoff, with $E^{(0)} = 0.1$.
2. Iteration: Updates are calculated using $\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} - \mathbf{J}_{\mathbf{F}}^{-1} \mathbf{F}(\mathbf{x}^{(k)})$.
3. Convergence: The process continues until the residual tolerance $\Vert{}\mathbf{F}\Vert{}_2 < 10^{-8}$ is reached.

### Reconciled Model Parameters

The following parameters ensure the numerical stability and physical accuracy of the integrated system:

Parameter	|$Symbol$	|$Value$	|$Units$	|$Physical / Numerical Role$
---|---|---|---|---
Regularization Floor	|$\delta$	|$10^{-6}$	|$\text{m}^2\,\text{s}^{-2}$	|Velocity scale smoothing constant
Buoyant Sink Efficiency	|$c_b$	|$0.50$	|—	|Dimensionless saturation coefficient
Shear Damping Coeff.	|$c_1$	|$1.80$	|$\text{s}\cdot\text{m}^{-2}$	|Turbulent feedback on mean shear
Timestep Safety Factor	|$\eta$	|$0.50$	|—	|RK4 step safety ($\Delta t_{\text{eff}} = \min(\Delta t, \frac{\eta \epsilon}{\vert{}\lambda_f\vert{}})$)

### Conclusion

The alignment of $f(E, S)$ and its derivative $f_E(E, S)$ ensures the system is mathematically closed. This synthesis allows for the reliable application of stiff ODE integration and precise tracking of critical fold transitions within the turbulent field.
