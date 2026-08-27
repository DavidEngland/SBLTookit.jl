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