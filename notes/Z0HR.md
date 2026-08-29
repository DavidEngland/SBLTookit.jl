The **Zero-Offset Hyperbolic Regularization (Z0HR)** scheme uses a brilliant mathematical trick to replace conditional piecewise logic (like `IF-ELSE` branches or `max`/`min` functions) with a single, branch-free, continuously differentiable (\(C^1\) or \(C^\infty\)) algebraic formulation.

It achieves this by decomposing the Richardson number (\(Ri\)) into smooth positive (\(Ri_+\)) and negative (\(Ri_-\)) components, and then **multiplying** the stable and unstable profile formulations together into a single unified equation.

---

### 1. How the Hyperbolic Coordinates (\(Ri_\pm\)) Work
Instead of using a sharp switch at \(Ri = 0\), Z0HR decomposes the Richardson number using smooth, hyperbolic coordinate transformations:
\[Ri_+ = \frac{Ri + \sqrt{Ri^2 + \epsilon^2}}{2} \tag{smooth max}\]
\[Ri_- = \frac{Ri - \sqrt{Ri^2 + \epsilon^2}}{2} \tag{smooth min}\]

The parameter \(\epsilon\) (typically set to \(10^{-3}\)) controls the smoothing window. The coordinate behavior shifts automatically depending on the sign of \(Ri\):
*   **When Stratification is Stable (\(Ri \gg 0\)):** \(Ri_+\) converges smoothly to the actual Richardson number (\(Ri\)), while \(Ri_-\) decays asymptotically to \(0\).
*   **When Stratification is Unstable (\(Ri \ll 0\)):** \(Ri_-\) converges smoothly to \(Ri\), while \(Ri_+\) decays asymptotically to \(0\).
*   **Near Neutrality (\(Ri \approx 0\)):** The coordinates smoothly transition and blend within an \(O(\epsilon)\) neighborhood, eliminating the discontinuous derivative "kink" at \(Ri = 0\).

---

### 2. The Unified Multiplication Strategy
The Z0HR stability functions for momentum (\(S_m\)) and heat (\(S_h\)) are structured as a **multiplication of a stable profile factor and an unstable profile factor**:

\[S_m(Ri) = \underbrace{\left[ f_{\text{stable}}(Ri) \right]^2}_{\text{Stable Profile}} \times \underbrace{\left(1 - B_{u,m} Ri_-\right)^{1/2}}_{\text{Unstable Profile}} \tag{1}\]

\[S_h(Ri) = \frac{1}{\alpha_\theta} \underbrace{\left[ f_{\text{stable}}(Ri) \right]^2}_{\text{Stable Profile}} \times \underbrace{\left(1 - B_{u,h} Ri_-\right)^{0.75}}_{\text{Unstable Profile}} \tag{2}\]

where the stable factor is regularized as \(f_{\text{stable}}(Ri) = \text{smooth\_max}\left(1 - \frac{Ri_+}{Ri_c}; \epsilon\right)\).

Because these two profiles are multiplied together, **the coordinate decomposition automatically activates and deactivates the correct regime** without requiring any logical branching.

---

### 3. Regime 1: Unstable Stratification (\(Ri < 0\))
When the atmosphere is unstable (convective conditions), we want the stability function to follow the classical power-law convective scaling. Here is how the multiplication achieves this:

1.  **Stable Term Deactivates:** Since \(Ri \ll 0\), the positive coordinate \(Ri_+ \to 0\). This means the stable factor simplifies to:
    \[f_{\text{stable}}(Ri) \to \text{smooth\_max}(1 - 0; \epsilon) \approx 1.0\]
    When squared, the stable profile term becomes a passive multiplier of \(1.0\).
2.  **Unstable Term Activates:** Meanwhile, \(Ri_- \to Ri\). The unstable profile factor becomes:
    \[\left(1 - B_u Ri\right)^p \quad (\text{where } p = 0.5 \text{ or } 0.75)\]
3.  **Result:** The multiplication yields the exact, asymptotically unbiased convective limit:
    \[S_m(Ri) \approx 1.0 \times \left(1 - B_{u,m} Ri\right)^{1/2}\]

---

### 4. Regime 2: Stable Stratification (\(Ri > 0\))
When the atmosphere is stably stratified, we want turbulent mixing to smoothly decay and shut off entirely once the critical Richardson number (\(Ri_c\)) is breached.

1.  **Unstable Term Deactivates:** Since \(Ri > 0\), the negative coordinate \(Ri_- \to 0\). This simplifies the unstable profile term to:
    \[\left(1 - B_u(0)\right)^p = 1.0^p = 1.0\]
    The unstable term is now a passive multiplier of \(1.0\).
2.  **Stable Term Activates:** Meanwhile, \(Ri_+ \to Ri\). The stable profile factor handles the decay:
    \[f_{\text{stable}}(Ri) \approx \text{smooth\_max}\left(1 - \frac{Ri}{Ri_c}; \epsilon\right)\]
3.  **Result:** The multiplication yields the exact stable cutoff:
    \[S_m(Ri) \approx \left[ \max\left(0, 1 - \frac{Ri}{Ri_c}\right) \right]^2 \times 1.0\]
    As \(Ri \to Ri_c\), the function smoothly approaches exact zero. For strong stability (\(Ri \ge Ri_c\)), it remains identically zero, successfully preventing unphysical spurious mixing in highly stratified nighttime layers.

---

### 5. Why This is Numerically Foolproof: Natural Positivity
In traditional codes, evaluating fractional powers like \((\cdot)^{1/2}\) or \((\cdot)^{0.75}\) on negative numbers during intermediate solver steps can crash the simulation with imaginary roots.

By using \(Ri_-\), Z0HR provides a **structural guarantee of mathematical positivity**:
*   Because \(Ri_-\) is strictly negative (\(Ri_- < 0\)) for all real \(Ri\), and the empirical constants \(B_u\) are strictly positive (\(B_u > 0\)), the term \(-B_u Ri_-\) is **always positive**.
*   Therefore, the radicand is structurally bounded:
    \[1 - B_u Ri_- \ge 1.0 \quad \forall Ri \in \mathbb{R}\]
*   This mathematical property guarantees that the fractional roots are always real, positive, and \(C^\infty\) smooth, completely eliminating the need for defensive, execution-slowing domain clamping.

---

