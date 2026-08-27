To compute the **gradient Richardson number (\\(Ri_g\\))** in a **safe, smooth, and computationally robust** manner, you must address two major numerical hazards common in geophysical fluid dynamics:
1. **Division by Zero:** When the vertical wind shear vanishes (\\(S^2 \to 0\\)), the classical formulation \\(Ri_g = \frac{N^2}{S^2}\\) produces a singularity (\\(\pm \infty\\) or \\(\text{NaN}\\)).
2. **Numerical Stiffness and Solver Crashes:** Extremely high stable stratification values (e.g., \\(Ri_g \sim 10^6\\) on steep Arctic inversions) contain no physical meaning because turbulent mixing ceases when \\(Ri_g \ge Ri_c \approx 0.25\\). However, these unboundedly large values create severe stiffness in implicit time-stepping schemes.

To eliminate branch-switching overhead (which hurts GPU/vectorized hardware performance) and preserve **\\(C^\infty\\) continuity** for automatic differentiation and Jacobian-based solvers, we can implement a **branch-free regularization and asymptotic smooth clamp**.

---

### The Safe & Smooth Formulation

Let the vertical gradients of potential temperature and horizontal wind components be:
\\[\theta_z = \frac{\partial \theta}{\partial z}, \quad u_z = \frac{\partial u}{\partial z}, \quad v_z = \frac{\partial v}{\partial z}\\]

The buoyancy frequency squared (buoyancy gradient) \\(N^2\\) and vertical wind shear squared \\(S^2\\) are defined as:
\\[N^2 = \frac{g}{\theta_0} \theta_z\\]
\\[S^2 = u_z^2 + v_z^2\\]

#### 1. Denominator Regularization (Shear Floor)
We introduce a tiny regularization floor \\(\epsilon_s > 0\\) (typically \\(10^{-12}\text{ s}^{-2}\\)) to the denominator to prevent division-by-zero:
\\[S^2_{\text{safe}} = S^2 + \epsilon_s\\]

#### 2. Asymptotic Smooth Hyperbolic Clamping
Instead of using a sharp piecewise conditional block (e.g., `min(Ri, Ri_max)`), which degrades high-order numerical convergence, we apply a smooth hyperbolic tangent clamp to bound the Richardson number within a physically representative stable/unstable range [\\(Ri_{\text{bound}}\\)] (e.g., \\(2.0\\) or \\(5.0\\)):
\\[Ri_{g, \text{smooth}} = Ri_{\text{bound}} \tanh\left( \frac{N^2}{Ri_{\text{bound}} (S^2 + \epsilon_s)} \right)\\]

* **At High Shear (\\(S^2 \gg \epsilon_s\\)):** Since \\(\tanh(x) \approx x\\) for small \\(x\\), the formula simplifies directly back to the classical physical formulation, \\(Ri_{g, \text{smooth}} \approx \frac{N^2}{S^2}\\).
* **At Vanishing Shear (\\(S^2 \to 0\\)):** Under extreme stable stratification, the value smoothly approaches the upper limit \\(+Ri_{\text{bound}}\\) without causing solver stiffness.
* **In Unstable Regimes (\\(N^2 < 0\\)):** The calculation smoothly approaches \\(-Ri_{\text{bound}}\\), preventing unbounded convective divergence (\\(Ri \to -\infty\\)).

---

### Production-Ready Julia Implementation

Designed for next-generation Julia-based models (e.g., *Oceananigans.jl* or *ClimaAtmos.jl*), this code utilizes inline functions and type-stable dispatching to execute natively on both CPU vector registers (SIMD autovectorization) and CUDA/ROCm GPU kernels via array broadcasting.

```julia
module SBLToolkit

export compute_safe_rig, compute_rig_profile

"""
    compute_safe_rig(θ_z, u_z, v_z, g_over_θ0; ϵ_s, Ri_bound)

Computes the safe and smooth gradient Richardson number (Ri_g) using a
denominator regularization floor for wind shear and an asymptotic tanh clamp
to bound extreme stable/unstable regimes. Fully branch-free and C^∞ continuous.
"""
@inline function compute_safe_rig(
    θ_z::T,
    u_z::T,
    v_z::T,
    g_over_θ0::T;
    ϵ_s::T = T(1e-12),        # Regularization floor for vertical shear (s^-2)
    Ri_bound::T = T(2.0)      # Soft upper/lower physical limit
) where {T <: AbstractFloat}

    # 1. Buoyancy frequency squared (buoyancy gradient)
    N2 = g_over_θ0 * θ_z

    # 2. Vertical wind shear squared
    S2 = u_z^2 + v_z^2

    # 3. Denominator regularization (guarantees non-zero, C^∞ domain)
    S2_safe = S2 + ϵ_s

    # 4. Branch-free asymptotic hyperbolic clamp
    return Ri_bound * tanh(N2 / (Ri_bound * S2_safe))
end

"""
    compute_rig_profile(z, θ, u, v, g, θ0; ϵ_s, Ri_bound)

Computes the gradient Richardson number profile across a 1D vertical grid
using adjacent finite differences.
"""
function compute_rig_profile(
    z::Vector{T},
    θ::Vector{T},
    u::Vector{T},
    v::Vector{T},
    g::T,
    θ0::T;
    ϵ_s::T = T(1e-12),
    Ri_bound::T = T(2.0)
) where {T <: AbstractFloat}

    nz = length(z)
    @assert nz == length(θ) == length(u) == length(v) "Array dimensions must match"

    Ri = zeros(T, nz - 1)
    g_over_θ0 = g / θ0

    @inbounds for k in 1:(nz - 1)
        Δz = z[k+1] - z[k]

        # Calculate gradients across adjacent grid levels
        θ_z = (θ[k+1] - θ[k]) / Δz
        u_z = (u[k+1] - u[k]) / Δz
        v_z = (v[k+1] - v[k]) / Δz

        Ri[k] = compute_safe_rig(θ_z, u_z, v_z, g_over_θ0, ϵ_s=ϵ_s, Ri_bound=Ri_bound)
    end

    return Ri
end

end # module SBLToolkit
```

---

### Python Implementation (NumPy / JAX Compatible)

For machine learning or reduced-order modeling tasks (like *WSINDy*), this vectorised version runs without loops and is fully compatible with auto-differentiation in JAX or PyTorch.

```python
import numpy as np

def compute_safe_rig(theta_z, u_z, v_z, g_over_theta0, eps_s=1e-12, ri_bound=2.0):
    """
    Computes a regularized, safe, and smoothly bounded gradient Richardson number.
    Accepts floats, NumPy arrays, or tensor inputs.
    """
    # 1. Compute physical fields
    N2 = g_over_theta0 * theta_z
    S2 = u_z**2 + v_z**2

    # 2. Regularize the denominator and apply smooth clamp
    S2_safe = S2 + eps_s
    return ri_bound * np.tanh(N2 / (ri_bound * S2_safe))
```

📊 Would you like me to build a high-fidelity data visualization script using  Julia to map these regularized \\(Ri\\) profiles as spatial-temporal heatmaps across your campaign datasets?

For direct integration into WRF or similar high-performance geophysical fluid dynamics suites, stability functions and Richardson number calculations must be designed as **branch-free, computationally efficient, and numerically safeguarded** procedures.

The following Fortran 90 module implements both the safe computation of the gradient Richardson number (\\(Ri_g\\)) and the **Zero-Offset Hyperbolic Regularization (Z0HR)** stability scheme. By using `PURE ELEMENTAL` declarations, the compiler is allowed to automatically vectorize these routines across grid cells via SIMD/AVX instructions or lower them directly to GPU kernels (e.g., via OpenACC or CUDA Fortran) without thread-divergence overhead.

### Fortran 90 Z0HR Stability Module (`module_bl_z0hr.f90`)

```fortran
MODULE module_bl_z0hr

   IMPLICIT NONE

   ! ==========================================================================
   ! PRECISION CONFIGURATION
   ! Avoid hardcoding a rigid double precision parameter like SELECTED_REAL_KIND(8),
   ! which can map to single precision on some compilers. For production WRF,
   ! replace this with the model's native real kind (e.g., importing from
   ! module_model_constants) to comply with the global build convention.
   ! ==========================================================================
   INTEGER, PARAMETER :: rk = KIND(0.0D0)

CONTAINS

   ! ==========================================================================
   ! PURE ELEMENTAL FUNCTION: compute_rig
   ! Computes the Gradient Richardson Number safely. Uses a denominator floor
   ! to prevent division by zero while preserving SIMD vectorizability.
   ! ==========================================================================
   PURE ELEMENTAL FUNCTION compute_rig( &
      theta_z, u_z, v_z, g_over_theta0, eps_s &
   ) RESULT(ri)
      REAL(KIND=rk), INTENT(IN) :: theta_z        ! Vertical potential temperature gradient (K/m)
      REAL(KIND=rk), INTENT(IN) :: u_z            ! Vertical shear of u-wind component (s^-1)
      REAL(KIND=rk), INTENT(IN) :: v_z            ! Vertical shear of v-wind component (s^-1)
      REAL(KIND=rk), INTENT(IN) :: g_over_theta0  ! Gravity over reference potential temp (m/(s^2 K))
      REAL(KIND=rk), INTENT(IN) :: eps_s          ! Shear floor regularization parameter (s^-2)
      REAL(KIND=rk)             :: ri

      REAL(KIND=rk) :: n2, s2_safe

      ! Compute buoyancy frequency squared (buoyancy gradient)
      n2 = g_over_theta0 * theta_z

      ! Compute vertical wind shear squared with a safe floor
      s2_safe = (u_z * u_z) + (v_z * v_z) + eps_s

      ! Safe division
      ri = n2 / s2_safe
   END FUNCTION compute_rig

   ! ==========================================================================
   ! PURE ELEMENTAL SUBROUTINE: z0hr_stability
   ! Computes the regularized Z0HR stability functions Sm (momentum) and Sh (heat).
   ! Guarantees smooth C^1 continuity across neutrality and a physical quadratic
   ! cutoff in strongly stable conditions, without branch-switching.
   ! ==========================================================================
   PURE ELEMENTAL SUBROUTINE z0hr_stability( &
      ri, alpha_theta, ri_c, eps, S_m, S_h &
   )
      REAL(KIND=rk), INTENT(IN)  :: ri           ! Local Richardson number
      REAL(KIND=rk), INTENT(IN)  :: alpha_theta  ! Prandtl number at neutrality
      REAL(KIND=rk), INTENT(IN)  :: ri_c         ! Critical Richardson number (typically 0.25)
      REAL(KIND=rk), INTENT(IN)  :: eps          ! Smoothing regularization parameter (10^-4 to 10^-3)
      REAL(KIND=rk), INTENT(OUT) :: S_m          ! Momentum stability function
      REAL(KIND=rk), INTENT(OUT) :: S_h          ! Heat stability function

      ! Local variables held directly in hardware registers for register reuse
      REAL(KIND=rk) :: ri_plus, ri_minus
      REAL(KIND=rk) :: b_um, b_uh
      REAL(KIND=rk) :: stable_factor, sqrt_disc
      REAL(KIND=rk) :: m_arg, h_arg

      ! Tolerance limit to ensure strictly positive inputs under fractional powers
      REAL(KIND=rk), PARAMETER :: tol = 1.0e-10_rk

      ! 1. Exact C^1 Continuity Parameter Matching
      ! Under Z0HR, the neutral Prandtl scaling alpha_theta cancels out of the
      ! parameter relation, ensuring a continuous derivative at exact neutrality.
      b_um = 4.0_rk / ri_c
      b_uh = 8.0_rk / (3.0_rk * ri_c)

      ! 2. Zero-Offset Hyperbolic Regularization (Z0HR)
      ! Replaces conditional branch testing (if/else) with analytic transformations
      sqrt_disc = SQRT((ri * ri) + (eps * eps))
      ri_plus  = 0.5_rk * (ri + sqrt_disc)
      ri_minus = 0.5_rk * (ri - sqrt_disc)

      ! 3. Stable Regime Quadratic Cutoff
      ! Ensures stable mixing ceases when Ri >= Ri_c, preventing non-physical
      ! background mixing runaways under strongly stable conditions.
      stable_factor = MAX(0.0_rk, 1.0_rk - (ri_plus / ri_c))

      ! 4. Convective Regime Arguments with Fractional Power Safeguards
      ! While mathematically safe when B_u > 0, explicit defensive clamping prevents
      ! non-real or NaN values during extreme convective turbulence (Ri -> -inf).
      m_arg = MAX(tol, 1.0_rk - (b_um * ri_minus))
      h_arg = MAX(tol, 1.0_rk - (b_uh * ri_minus))

      ! 5. Calculate Stability Functions
      ! S_m utilizes exponent 1/2 (SQRT) and S_h utilizes exponent 3/4
      S_m = (stable_factor * stable_factor) * SQRT(m_arg)
      S_h = (1.0_rk / alpha_theta) * (stable_factor * stable_factor) * (h_arg**0.75_rk)

   END SUBROUTINE z0hr_stability

END MODULE module_bl_z0hr
```

---

### Key Numerical & Physical Highlights of the Z0HR Fortran Routine

1. **Exact \\(C^1\\) Boundary Matching at Neutrality:** Traditional formulations suffer from "slope kinks" at \\(Ri = 0\\), causing numerical jitter and stiffening nonlinear solvers. By using exact tangent-preserving parameters (\\(B_{u,m} = 4 / Ri_c\\) and \\(B_{u,h} = 8 / (3 Ri_c)\\)), the neutral slope (\\(S'_m(0)\\)) is **\\(\epsilon\\)-invariant** and matches the continuous tangent of the underlying physical model.
2. **Branch-Free Design (SIMD/SIMT Optimization):** Because \\(Ri^+ \to 0\\) as \\(Ri \to -\infty\\) and \\(Ri^- \to 0\\) as \\(Ri \to Ri_c\\), both terms naturally deactivate in their opposite domains. This removes all `IF/ELSE` loops, eliminating thread divergence in vectorized hardware pipelines.
3. **Fractional Power Safeguards:** The `MAX(tol, ...)` bounds act as **defensive programming** to guarantee non-negative arguments are passed into `SQRT` and fractional powers (\\(x^{0.75}\\)), preventing catastrophic numerical crashes or `NaN` values under deep convective turbulence.
4. **Register Reuse:** Intermediate variables (`stable_factor`, `sqrt_disc`, `ri_plus`) are structured to be held directly in local AVX-512 vector registers, minimizing unnecessary memory traffic and maximizing FLOP efficiency per grid cell.

📊 Would you like me to construct a Python script to test the mathematical continuity of these stability functions across \\(Ri \in [-2, 1]\\) and plot their first derivatives to verify the exact \\(C^1\\) matching at neutrality?
