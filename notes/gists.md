### Gist 1: Core Physics Module (`module_bl_cz0hr.F90`)

**Gist Description:** `C-Z0HR regularized Richardson number & C^inf smooth stability function module for WRF/MPAS PBL physics (module_bl_cz0hr.F90)`

```fortran
!===============================================================================
! MODULE: module_bl_cz0hr
!
! PURPOSE:
!   Provides Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR)
!   for Planetary Boundary Layer (PBL) schemes in NWP models (WRF / MPAS).
!
! TECHNICAL SUMMARY:
!   - Replaces C^0 piecewise clamps with C^inf hyperbolic softplus functions.
!   - Modifies input Gradient Richardson Number (Ri_g) via local curvature C.
!   - Attenuates threshold sensitivity (dRi_reg/dRi -> 1/(1+alpha) at Ri_c).
!   - Preserves unmodified legacy K(Ri) closures without retuning parameters.
!   - PURE ELEMENTAL routines ensure zero heap allocation, SIMD vectorization,
!     and compatibility with OpenACC / GPU offloading.
!
! KEYWORDS / SEARCH TAGS:
!   WRF, MPAS, PBL Physics, Richardson Number, Boundary Layer, Fortran 90,
!   Regularization, Zero-Offset Hyperbolic Regularization, Turbulence Closure
!===============================================================================
MODULE module_bl_cz0hr

  IMPLICIT NONE

  ! Precision control (Standard double precision for Jacobian stability)
  INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(12, 60)

CONTAINS

  !-----------------------------------------------------------------------------
  ! Smooth C^inf Hyperbolic Distance Metric
  !   Phi_eps(x) = SQRT(x^2 + eps^2)
  !-----------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION phi_eps(x, eps_hyper) RESULT(res)
    REAL(rk), INTENT(IN) :: x, eps_hyper
    REAL(rk)             :: res

    res = SQRT(x * x + eps_hyper * eps_hyper)
  END FUNCTION phi_eps

  !-----------------------------------------------------------------------------
  ! Smooth C^inf Hyperbolic Softplus Operator
  !   softplus(g) = 0.5 * (g + SQRT(g^2 + eps^2))
  !-----------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION smooth_softplus(g_val, eps_hyper) RESULT(res)
    REAL(rk), INTENT(IN) :: g_val, eps_hyper
    REAL(rk)             :: res

    res = 0.5_rk * (g_val + SQRT(g_val * g_val + eps_hyper * eps_hyper))
  END FUNCTION smooth_softplus

  !-----------------------------------------------------------------------------
  ! Core C-Z0HR Mapping Operator
  !   Transforms raw Ri_g into regularized Ri_reg based on profile curvature.
  !-----------------------------------------------------------------------------
  PURE ELEMENTAL SUBROUTINE cz0hr_transform_ri( &
      Ri_raw, d2Ri_dz2, dz, Ri_c, alpha, eps_c, eps_hyper, &
      Ri_reg, D_alpha, C_scale)

    ! Inputs
    REAL(rk), INTENT(IN)  :: Ri_raw     ! Unregularized Gradient Richardson Number
    REAL(rk), INTENT(IN)  :: d2Ri_dz2   ! 2nd spatial vertical derivative (d²Ri/dz²)
    REAL(rk), INTENT(IN)  :: dz         ! Vertical grid spacing [m]
    REAL(rk), INTENT(IN)  :: Ri_c       ! Critical Richardson threshold (e.g., 0.20 or 0.25)
    REAL(rk), INTENT(IN)  :: alpha      ! Damping control factor (default: 2.0)
    REAL(rk), INTENT(IN)  :: eps_c      ! Softening term scaling factor (default: 0.05)
    REAL(rk), INTENT(IN)  :: eps_hyper  ! Hyperbolic smoothing parameter (default: 1e-3)

    ! Outputs
    REAL(rk), INTENT(OUT) :: Ri_reg     ! Regularized Richardson number input for K(Ri)
    REAL(rk), INTENT(OUT) :: D_alpha    ! Partial slope (dRi_reg / dRi)
    REAL(rk), INTENT(OUT) :: C_scale    ! Evaluated spatial curvature scale C

    ! Local variables
    REAL(rk) :: x, phi, num_reg, den_reg, num_D, den_D

    ! 1. Distance from critical threshold
    x = Ri_raw - Ri_c

    ! 2. Hyperbolic distance metric
    phi = phi_eps(x, eps_hyper)

    ! 3. Discretized spatial curvature scale C
    C_scale = 0.5_rk * ABS(d2Ri_dz2) * (dz * dz) + (eps_c * Ri_c)

    ! 4. C-Z0HR Richardson mapping
    num_reg = phi + C_scale
    den_reg = phi + (1.0_rk + alpha) * C_scale
    Ri_reg  = Ri_c + x * (num_reg / den_reg)

    ! 5. Partial slope evaluation: D_alpha = d(Ri_reg)/d(Ri_raw)
    num_D   = (phi * phi) + 2.0_rk * (1.0_rk + alpha) * C_scale * phi + &
              (1.0_rk + alpha) * (C_scale * C_scale)
    den_D   = den_reg * den_reg
    D_alpha = num_D / den_D

  END SUBROUTINE cz0hr_transform_ri

  !-----------------------------------------------------------------------------
  ! C^inf Smooth Momentum Stability Function S_m(Ri)
  !-----------------------------------------------------------------------------
  PURE ELEMENTAL FUNCTION cz0hr_stability_sm(Ri_reg, Ri_c, eps_hyper) RESULT(Sm)
    REAL(rk), INTENT(IN) :: Ri_reg, Ri_c, eps_hyper
    REAL(rk)             :: Sm
    REAL(rk)             :: g_raw, soft_g

    g_raw  = 1.0_rk - (Ri_reg / Ri_c)
    soft_g = smooth_softplus(g_raw, eps_hyper)
    Sm     = soft_g * soft_g

  END FUNCTION cz0hr_stability_sm

END MODULE module_bl_cz0hr

```

---

### Gist 2: WRF 1D Column PBL Driver Integration (`cz0hr_column_driver.F90`)

**Gist Description:** `WRF-compatible 1D column physics driver implementing C-Z0HR for vertical diffusion (cz0hr_column_driver.F90)`

```fortran
!===============================================================================
! SUBROUTINE: cz0hr_column_driver
!
! EXAMPLE INTEGRATION FOR WRF PBL SCHEMES (e.g., YSU, MYNN, or Louis extensions)
!===============================================================================
SUBROUTINE cz0hr_column_driver( &
    kts, kte, dz, u, v, th, qv, &
    Ri_c, alpha_param, Km_out, Kh_out)

  USE module_bl_cz0hr, ONLY : rk, cz0hr_transform_ri, cz0hr_stability_sm
  IMPLICIT NONE

  ! Subroutine Arguments
  INTEGER,  INTENT(IN)  :: kts, kte       ! Vertical start/end indices
  REAL(rk), INTENT(IN)  :: dz             ! Uniform vertical grid spacing [m]
  REAL(rk), DIMENSION(kts:kte), INTENT(IN)  :: u, v, th, qv ! State variables
  REAL(rk), INTENT(IN)  :: Ri_c           ! Critical Richardson threshold
  REAL(rk), INTENT(IN)  :: alpha_param    ! Damping control parameter (e.g., 2.0)
  REAL(rk), DIMENSION(kts:kte), INTENT(OUT) :: Km_out, Kh_out ! Diffusivities

  ! Local Arrays
  REAL(rk), DIMENSION(kts:kte) :: S_shear, N2_buoy, Ri_raw, d2Ri_dz2
  REAL(rk), DIMENSION(kts:kte) :: Ri_reg, D_alpha, C_scale, l_m
  REAL(rk) :: du_dz, dv_dz, dth_dz, Sm_val, g_acc, th_0
  INTEGER  :: k

  PARAMETER (g_acc = 9.81_rk, th_0 = 290.0_rk)

  ! 1. Compute local vertical shear (S) and buoyancy frequency (N^2)
  DO k = kts, kte-1
    du_dz  = (u(k+1) - u(k)) / dz
    dv_dz  = (v(k+1) - v(k)) / dz
    dth_dz = (th(k+1) - th(k)) / dz

    S_shear(k) = SQRT(du_dz*du_dz + dv_dz*dv_dz + 1.0e-8_rk)
    N2_buoy(k) = (g_acc / th_0) * dth_dz
    Ri_raw(k)  = N2_buoy(k) / (S_shear(k) * S_shear(k))
  END DO
  Ri_raw(kte) = Ri_raw(kte-1)
  S_shear(kte) = S_shear(kte-1)

  ! 2. Compute 2nd spatial derivative of Richardson profile: d²Ri/dz²
  DO k = kts+1, kte-1
    d2Ri_dz2(k) = (Ri_raw(k+1) - 2.0_rk * Ri_raw(k) + Ri_raw(k-1)) / (dz * dz)
  END DO
  d2Ri_dz2(kts) = d2Ri_dz2(kts+1)
  d2Ri_dz2(kte) = d2Ri_dz2(kte-1)

  ! 3. Apply C-Z0HR Regularization Mapping across column
  DO k = kts, kte
    CALL cz0hr_transform_ri( &
        Ri_raw(k), d2Ri_dz2(k), dz, Ri_c, alpha_param, &
        0.05_rk, 1.0e-3_rk, Ri_reg(k), D_alpha(k), C_scale(k))

    ! Mix length scale (Blackadar formulation)
    l_m(k) = 15.0_rk * (REAL(k, rk)*dz / (REAL(k, rk)*dz + 20.0_rk))

    ! Evaluate C^inf smooth stability function
    Sm_val = cz0hr_stability_sm(Ri_reg(k), Ri_c, 1.0e-3_rk)

    ! Evaluate unmodified K_m closure
    Km_out(k) = (l_m(k) * l_m(k)) * S_shear(k) * Sm_val
    Kh_out(k) = Km_out(k) ! Pr_t = 1.0 assumption
  END DO

END SUBROUTINE cz0hr_column_driver

```

---

### GitHub Gist SEO & Discoverability Strategy

To maximize search visibility on GitHub and search engines (Google, Bing) when WRF/MPAS developers search for stability scheme fixes or Richardson number smoothers:

**Gist Titles**

* `WRF / MPAS Fortran 90 C-Z0HR Richardson Number Regularization Scheme`
* `Branch-Free C^inf Smooth PBL Stability Closure for WRF Vertical Diffusion`

**Tags & Topic Hashtags**
`#WRF` `#MPAS` `#Fortran90` `#PlanetaryBoundaryLayer` `#NumericalWeatherPrediction` `#RichardsonNumber` `#AtmosphericPhysics` `#OpenACC` `#SIMD`

**Key Search Query Keywords Embedded in Header Comments**

* WRF PBL physics module
* Richardson number cutoff singularity fix
* Branch-free smooth stability function
* Non-differentiable $C^0$ kink elimination
* Newton-Raphson convergence acceleration for implicit solvers
* OpenACC/OpenMP PURE ELEMENTAL Fortran 90 subroutine
