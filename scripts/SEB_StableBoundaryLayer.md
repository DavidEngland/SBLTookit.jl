---
title: SEB Stable Boundary Layer Model
author: Refactor of seb.vba to SEB_StableBoundaryLayer.bas
date: 2024-06-05
---

**What the model does**

It's a single-column nocturnal boundary-layer model. Picture a vertical column of air from 1.5 m up to 200 m, discretized into `Nz` levels:

1. **Surface energy balance**: each time step, the surface temperature `Ts` is updated implicitly from four terms — its own heat storage (`Cs`), a fixed net radiative cooling (`Rn0`, negative = net loss at night), turbulent exchange with the lowest air level, and conductive coupling to a deep soil reservoir at `T_deep`.
2. **Turbulence closure**: how efficiently heat mixes between the surface and the air (the eddy diffusivity `Kh`) depends on how *stable* the air is, measured by the Richardson number `Ri`. Under strong stability, turbulence should shut off — but exactly how it shuts off is the crux of this comparison:
   - **Standard closure**: `Kh` scales with `(1 − Ri/Ri_c)²`, clipped to zero once `Ri` exceeds the critical value `Ri_c`. This clipping has a kink (discontinuous derivative) at `Ri_c`, which is a known source of grid-dependent, sometimes runaway cooling in stable boundary-layer models.
   - **C-Z0HR closure**: a smoothed version of the same function that blends toward the raw Richardson number away from `Ri_c` and applies a curvature-aware smoothing near it, removing the hard cutoff.
3. **Implicit vertical diffusion**: each step, the whole column's temperature profile is relaxed via a tridiagonal (Thomas-algorithm) solve, coupled to the surface via `Kh(1)`.
4. **Benchmark**: `RunSEBBenchmark` runs both closures at four grid resolutions (25, 50, 100, 200 levels) and prints final surface temperature, total cooling, sensible heat flux, and near-surface `Kh` — the point being to see whether the standard scheme's results drift with grid resolution while C-Z0HR's stay stable.

**How to use it**

1. Import `SEB_StableBoundaryLayer.bas` into a VBA project (Excel/Access: Alt+F11 → File → Import File…).
2. Open the Immediate Window (Ctrl+G) so you can see `Debug.Print` output.
3. Run `RunSEBBenchmark` (press F5 with the cursor in that sub, or type `RunSEBBenchmark` in the Immediate Window and hit Enter).
4. To run a single case yourself instead of the full benchmark:

```vba
Dim cfg As SEBModelConfig
cfg = GetDefaultConfig(100)   ' 100 grid levels

Dim Ts As Double, dTs As Double, H0 As Double, Kh As Double
RunSEBSimulation True, cfg, Ts, dTs, H0, Kh   ' True = use C-Z0HR closure
Debug.Print "Final Ts="; Ts; " H0="; H0
```

1. To change the scenario (e.g. stronger radiative cooling, different mixing length), edit the fields on the `cfg` returned by `GetDefaultConfig`, or start from a fresh `SEBModelConfig` and set every field yourself before calling `RunSEBSimulation`.

---

```vba
Option Explicit

'===============================================================================
' MODULE: SEB_StableBoundaryLayer
'
' Single-column surface energy balance (SEB) model coupled to an implicit
' vertical diffusion scheme for potential temperature. Used to compare two
' eddy-diffusivity closures under stable (nocturnal) stratification:
'
'   - "Standard" : classic Richardson-number stability function, hard-clipped
'                  to zero once Ri exceeds the critical Richardson number Ri_c
'   - "C-Z0HR"   : a smooth ("C-infinity") regularization of the same closure
'                  that removes the hard cutoff at Ri_c, avoiding the
'                  numerical stiffness/discontinuity that produces grid-
'                  dependent results with the standard scheme
'
' RunSEBBenchmark runs both closures across several grid resolutions and
' prints a side-by-side comparison so the grid sensitivity of each scheme
' can be inspected directly.
'===============================================================================

' ---- Physical / numerical constants (not user-configurable) ----
Private Const G_GRAVITY As Double = 9.81   ' gravitational acceleration, m/s^2
Private Const THETA_REF As Double = 285#   ' reference potential temperature, K
Private Const Z_SURFACE As Double = 1.5    ' height of lowest grid level, m
Private Const THETA_INIT As Double = 285#  ' initial column temperature, K

Public Type SEBModelConfig
    ' --- Grid / time ---
    Nz As Long          ' number of vertical levels
    z_max As Double     ' top of domain, m
    dt As Double        ' time step, s
    t_end As Double     ' simulation length, s

    ' --- Surface energy balance ---
    Cs As Double        ' surface heat capacity, J/(m^2 K)
    Ks As Double        ' deep-soil coupling conductance, W/(m^2 K)
    T_deep As Double    ' deep-soil temperature, K
    Rn0 As Double       ' net radiation forcing, W/m^2
    rho_cp As Double    ' air density * specific heat, J/(m^3 K)

    ' --- Turbulence closure ---
    Ri_c As Double      ' critical Richardson number
    alpha As Double     ' C-Z0HR regularization sharpness parameter
    epsilon_c As Double ' C-Z0HR floor for the curvature scale
    eps_hyper As Double ' smoothing parameter for Phi_eps(), the smoothed |x|
    l_m As Double       ' mixing length, m
    S_eff As Double     ' effective shear, 1/s
End Type

' ------------------------------------------------------------------
' Returns a config populated with the model's default scenario values.
' ------------------------------------------------------------------
Public Function GetDefaultConfig(Optional ByVal grid_levels As Long = 50) As SEBModelConfig
    Dim cfg As SEBModelConfig

    cfg.Nz = grid_levels
    cfg.z_max = 200#
    cfg.dt = 60#
    cfg.t_end = 12# * 3600#

    cfg.Cs = 200000#
    cfg.Ks = 2#
    cfg.T_deep = 285#
    cfg.Rn0 = -50#
    cfg.rho_cp = 1200#

    cfg.Ri_c = 0.2
    cfg.alpha = 2#
    cfg.epsilon_c = 0.05
    cfg.eps_hyper = 0.001
    cfg.l_m = 15#
    cfg.S_eff = 0.015

    GetDefaultConfig = cfg
End Function

' ------------------------------------------------------------------
' Smoothed absolute value: sqrt(x^2 + eps^2). As eps_h -> 0 this
' approaches |x|, but stays differentiable everywhere.
' ------------------------------------------------------------------
Private Function Phi_eps(ByVal x As Double, ByVal eps_h As Double) As Double
    Phi_eps = Sqr(x * x + eps_h * eps_h)
End Function

' ------------------------------------------------------------------
' Builds the vertical grid: Nz levels evenly spaced from Z_SURFACE to
' cfg.z_max. Returns the level heights and the (uniform) spacing dz.
' ------------------------------------------------------------------
Private Function BuildHeightGrid(ByVal Nz As Long, ByVal z_max As Double, ByRef dz As Double) As Double()
    Dim z() As Double
    ReDim z(1 To Nz)
    Dim i As Long

    dz = (z_max - Z_SURFACE) / (Nz - 1)
    For i = 1 To Nz
        z(i) = Z_SURFACE + (i - 1) * dz
    Next i

    BuildHeightGrid = z
End Function

' ------------------------------------------------------------------
' Vertical gradient of potential temperature at each level, using the
' surface value Ts for the first level and simple differences above.
' ------------------------------------------------------------------
Private Function ComputeThetaGradient(ByRef theta() As Double, ByVal Ts As Double, _
                                       ByRef z() As Double, ByVal dz As Double, ByVal Nz As Long) As Double()
    Dim dtheta_dz() As Double
    ReDim dtheta_dz(1 To Nz)
    Dim i As Long

    dtheta_dz(1) = (theta(1) - Ts) / z(1)
    For i = 2 To Nz
        dtheta_dz(i) = (theta(i) - theta(i - 1)) / dz
    Next i

    ComputeThetaGradient = dtheta_dz
End Function

' ------------------------------------------------------------------
' Gradient Richardson number at each level (unregularized/raw value).
' ------------------------------------------------------------------
Private Function ComputeRichardson(ByRef dtheta_dz() As Double, ByVal S_eff As Double, ByVal Nz As Long) As Double()
    Dim Ri_raw() As Double
    ReDim Ri_raw(1 To Nz)
    Dim i As Long

    For i = 1 To Nz
        Ri_raw(i) = (G_GRAVITY / THETA_REF) * dtheta_dz(i) / (S_eff ^ 2)
    Next i

    ComputeRichardson = Ri_raw
End Function

' ------------------------------------------------------------------
' Second derivative of the raw Richardson profile (central differences,
' one-sided/copied at the boundaries). Used only by the C-Z0HR closure
' to estimate a local curvature scale for its regularization.
' ------------------------------------------------------------------
Private Function ComputeRichardsonCurvature(ByRef Ri_raw() As Double, ByVal dz As Double, ByVal Nz As Long) As Double()
    Dim d2Ri_dz2() As Double
    ReDim d2Ri_dz2(1 To Nz)
    Dim i As Long

    d2Ri_dz2(1) = 0#
    For i = 2 To Nz - 1
        d2Ri_dz2(i) = (Ri_raw(i + 1) - 2# * Ri_raw(i) + Ri_raw(i - 1)) / (dz ^ 2)
    Next i
    d2Ri_dz2(Nz) = d2Ri_dz2(Nz - 1)

    ComputeRichardsonCurvature = d2Ri_dz2
End Function

' ------------------------------------------------------------------
' Standard (hard-clipped) Louis-type stability function:
'   Kh = l_m^2 * S_eff * max(1 - Ri/Ri_c, 0)^2
' Discontinuous derivative at Ri = Ri_c is the source of this scheme's
' grid sensitivity.
' ------------------------------------------------------------------
Private Function ComputeDiffusivity_Standard(ByRef Ri_raw() As Double, ByRef cfg As SEBModelConfig, ByVal Nz As Long) As Double()
    Dim Kh() As Double
    ReDim Kh(1 To Nz)
    Dim i As Long, g_val As Double

    For i = 1 To Nz
        g_val = 1# - Ri_raw(i) / cfg.Ri_c
        If g_val < 0# Then g_val = 0#
        Kh(i) = (cfg.l_m ^ 2) * cfg.S_eff * (g_val ^ 2)
    Next i

    ComputeDiffusivity_Standard = Kh
End Function

' ------------------------------------------------------------------
' C-Z0HR closure: a smooth (C-infinity) regularization of the same
' stability function that removes the hard cutoff at Ri_c. Locally
' blends toward the raw Richardson number away from Ri_c, and toward a
' smoothed-cutoff behavior near Ri_c, scaled by the profile's local
' curvature (d2Ri_dz2) so the regularization tightens where the raw
' profile is smooth and loosens where it is not.
' ------------------------------------------------------------------
Private Function ComputeDiffusivity_CZ0HR(ByRef Ri_raw() As Double, ByRef cfg As SEBModelConfig, ByVal dz As Double, ByVal Nz As Long) As Double()
    Dim Kh() As Double
    ReDim Kh(1 To Nz)
    Dim d2Ri_dz2() As Double
    Dim i As Long
    Dim x_val As Double, Phi As Double, C_scale As Double
    Dim Ri_reg As Double, g_reg As Double, softplus_val As Double

    d2Ri_dz2 = ComputeRichardsonCurvature(Ri_raw, dz, Nz)

    For i = 1 To Nz
        x_val = Ri_raw(i) - cfg.Ri_c
        Phi = Phi_eps(x_val, cfg.eps_hyper)
        C_scale = 0.5 * Abs(d2Ri_dz2(i)) * (dz ^ 2) + cfg.epsilon_c * cfg.Ri_c

        Ri_reg = cfg.Ri_c + x_val * ((Phi + C_scale) / (Phi + (1# + cfg.alpha) * C_scale))

        g_reg = 1# - Ri_reg / cfg.Ri_c
        softplus_val = 0.5 * (g_reg + Phi_eps(g_reg, cfg.eps_hyper))
        Kh(i) = (cfg.l_m ^ 2) * cfg.S_eff * (softplus_val ^ 2)
    Next i

    ComputeDiffusivity_CZ0HR = Kh
End Function

' ------------------------------------------------------------------
' Dispatches to the requested closure and returns Kh(1..Nz).
' ------------------------------------------------------------------
Private Function ComputeDiffusivity(ByVal use_cz0hr As Boolean, ByRef Ri_raw() As Double, _
                                     ByRef cfg As SEBModelConfig, ByVal dz As Double, ByVal Nz As Long) As Double()
    If use_cz0hr Then
        ComputeDiffusivity = ComputeDiffusivity_CZ0HR(Ri_raw, cfg, dz, Nz)
    Else
        ComputeDiffusivity = ComputeDiffusivity_Standard(Ri_raw, cfg, Nz)
    End If
End Function

' ------------------------------------------------------------------
' Implicit (backward-Euler) update of the surface temperature from the
' surface energy balance: heat storage + net radiation + turbulent
' exchange with level 1 + conductive exchange with the deep soil.
' Also returns the resulting sensible heat flux H0 (positive upward;
' H0 < 0 indicates downward/into-the-surface heat exchange).
' ------------------------------------------------------------------
Private Sub UpdateSurfaceTemperature(ByRef Ts As Double, ByVal theta1 As Double, ByVal Kh1 As Double, _
                                      ByVal z1 As Double, ByRef cfg As SEBModelConfig, ByRef H0 As Double)
    Dim gamma_coupling As Double, beta_denom As Double

    gamma_coupling = cfg.rho_cp * Kh1 / z1
    beta_denom = (cfg.Cs / cfg.dt) + gamma_coupling + cfg.Ks

    Ts = ((cfg.Cs / cfg.dt) * Ts + cfg.Rn0 + gamma_coupling * theta1 + cfg.Ks * cfg.T_deep) / beta_denom
    H0 = cfg.rho_cp * Kh1 * (Ts - theta1) / z1
End Sub

' ------------------------------------------------------------------
' Thomas algorithm (tridiagonal solve). a = sub-diagonal, b = main
' diagonal, c = super-diagonal, d = right-hand side; result in x.
' ------------------------------------------------------------------
Private Sub SolveTridiagonal(ByRef a() As Double, ByRef b() As Double, ByRef c() As Double, _
                             ByRef d() As Double, ByRef x() As Double, ByVal N As Long)
    Dim i As Long
    Dim c_prime() As Double, d_prime() As Double
    ReDim c_prime(1 To N) As Double
    ReDim d_prime(1 To N) As Double
    Dim m As Double

    c_prime(1) = c(1) / b(1)
    d_prime(1) = d(1) / b(1)

    For i = 2 To N
        m = b(i) - a(i) * c_prime(i - 1)
        If i < N Then c_prime(i) = c(i) / m
        d_prime(i) = (d(i) - a(i) * d_prime(i - 1)) / m
    Next i

    x(N) = d_prime(N)
    For i = N - 1 To 1 Step -1
        x(i) = d_prime(i) - c_prime(i) * x(i + 1)
    Next i
End Sub

' ------------------------------------------------------------------
' Assembles the tridiagonal system for one implicit diffusion step of
' the column and solves it in place (theta is updated on return).
' Level 1 is coupled to the surface via r_s; interior levels use the
' standard implicit-diffusion stencil; the top level has a zero-flux
' (Neumann) upper boundary.
' ------------------------------------------------------------------
Private Sub DiffuseColumnImplicit(ByRef theta() As Double, ByRef Kh() As Double, ByVal Ts As Double, _
                                   ByVal z1 As Double, ByVal dz As Double, ByVal dt As Double, ByVal Nz As Long)
    Dim sub_diag() As Double, main_diag() As Double, sup_diag() As Double, rhs() As Double
    ReDim sub_diag(1 To Nz): ReDim main_diag(1 To Nz): ReDim sup_diag(1 To Nz): ReDim rhs(1 To Nz)
    Dim i As Long, r_i As Double

    Dim r_s As Double: r_s = Kh(1) * dt / (z1 ^ 2)
    Dim r_v As Double: r_v = Kh(1) * dt / (dz ^ 2)

    main_diag(1) = 1# + r_s + r_v
    sup_diag(1) = -r_v
    sub_diag(1) = 0#
    rhs(1) = theta(1) + r_s * Ts

    For i = 2 To Nz - 1
        r_i = Kh(i) * dt / (dz ^ 2)
        sub_diag(i) = -r_i
        main_diag(i) = 1# + 2# * r_i
        sup_diag(i) = -r_i
        rhs(i) = theta(i)
    Next i

    Dim r_N As Double: r_N = Kh(Nz) * dt / (dz ^ 2)
    sub_diag(Nz) = -r_N
    main_diag(Nz) = 1# + r_N
    sup_diag(Nz) = 0#
    rhs(Nz) = theta(Nz)

    SolveTridiagonal sub_diag, main_diag, sup_diag, rhs, theta, Nz
End Sub

' ------------------------------------------------------------------
' Runs the full time-marching SEB + diffusion simulation for one
' closure choice and returns the final diagnostics.
' ------------------------------------------------------------------
Public Sub RunSEBSimulation(ByVal use_cz0hr As Boolean, ByRef cfg As SEBModelConfig, _
                            ByRef final_Ts As Double, ByRef delta_Ts As Double, _
                            ByRef final_H0 As Double, ByRef final_Kh As Double)
    Dim Nz As Long: Nz = cfg.Nz
    Dim dz As Double
    Dim z() As Double: z = BuildHeightGrid(Nz, cfg.z_max, dz)

    Dim theta() As Double: ReDim theta(1 To Nz)
    Dim i As Long
    For i = 1 To Nz
        theta(i) = THETA_INIT
    Next i
    Dim Ts As Double: Ts = THETA_INIT

    Dim Nt As Long: Nt = CLng(cfg.t_end / cfg.dt)
    Dim n As Long
    Dim H0 As Double
    Dim dtheta_dz() As Double, Ri_raw() As Double, Kh() As Double

    For n = 1 To Nt
        dtheta_dz = ComputeThetaGradient(theta, Ts, z, dz, Nz)
        Ri_raw = ComputeRichardson(dtheta_dz, cfg.S_eff, Nz)
        Kh = ComputeDiffusivity(use_cz0hr, Ri_raw, cfg, dz, Nz)

        UpdateSurfaceTemperature Ts, theta(1), Kh(1), z(1), cfg, H0

        DiffuseColumnImplicit theta, Kh, Ts, z(1), dz, cfg.dt, Nz
    Next n

    final_Ts = Ts
    delta_Ts = Ts - THETA_INIT
    final_H0 = H0
    final_Kh = Kh(1)
End Sub

' ------------------------------------------------------------------
' Prints one formatted results row for the benchmark table.
' ------------------------------------------------------------------
Private Sub PrintSchemeRow(ByVal label As String, ByVal Ts As Double, ByVal dTs As Double, _
                            ByVal H0 As Double, ByVal Kh As Double)
    Debug.Print Format(label, "@@@@@@@@@@@@@@@") & " | " & _
                Format(Format(Ts, "0.00"), "@@@@@@@@@@@@") & " | " & _
                Format(Format(dTs, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                Format(Format(H0, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                Format(Format(Kh, "0.00E+00"), "@@@@@@@@@@@@@@@")
End Sub

' ------------------------------------------------------------------
' Runs both closures across a range of grid resolutions and prints a
' side-by-side comparison, so grid sensitivity can be inspected.
' ------------------------------------------------------------------
Public Sub RunSEBBenchmark()
    Dim grid_sizes As Variant
    grid_sizes = Array(25, 50, 100, 200)

    Dim idx As Long
    For idx = LBound(grid_sizes) To UBound(grid_sizes)
        Dim Nz As Long: Nz = grid_sizes(idx)
        Dim cfg As SEBModelConfig: cfg = GetDefaultConfig(Nz)

        Dim Ts_u As Double, dTs_u As Double, H0_u As Double, Kh_u As Double
        Dim Ts_c As Double, dTs_c As Double, H0_c As Double, Kh_c As Double

        RunSEBSimulation False, cfg, Ts_u, dTs_u, H0_u, Kh_u
        RunSEBSimulation True, cfg, Ts_c, dTs_c, H0_c, Kh_c

        Debug.Print "========================================================================="
        Debug.Print " GRID RESOLUTION AUDIT: Nz = " & Nz & " (dz = " & Format((200 - 1.5) / (Nz - 1), "0.00") & " m)"
        Debug.Print "========================================================================="
        Debug.Print "Scheme         | Final Ts (K) | Cooling dTs (K) | Final H0 (W/m2) | Final Kh (m2/s)"
        Debug.Print "-------------------------------------------------------------------------"
        PrintSchemeRow "Unregularized", Ts_u, dTs_u, H0_u, Kh_u
        PrintSchemeRow "C-Z0HR Reg.", Ts_c, dTs_c, H0_c, Kh_c
        Debug.Print "=========================================================================" & vbCrLf
    Next idx
End Sub

```

---

Refactor of seb.vba to SEB_StableBoundaryLayer.bas

**What changed in the refactor**

The math is byte-for-byte identical to your original — I didn't touch a single formula. What I did was break the two big procedures (`RunSEBSimulation`'s monolithic loop and `RunSEBBenchmark`'s inline printing) into named, single-purpose functions:

- `BuildHeightGrid`, `ComputeThetaGradient`, `ComputeRichardson`, `ComputeRichardsonCurvature` — pull the per-step array computations out into functions that return the array they build
- `ComputeDiffusivity_Standard` / `ComputeDiffusivity_CZ0HR` / `ComputeDiffusivity` (dispatcher) — the two closures are now separately readable and testable instead of living inside one big `If/Else`
- `UpdateSurfaceTemperature`, `DiffuseColumnImplicit` — the SEB update and the tridiagonal-system assembly are isolated from the closure logic
- `PrintSchemeRow` — removes the duplicated `Format(...)` block that appeared twice in the benchmark
- Named constants (`G_GRAVITY`, `THETA_REF`, `Z_SURFACE`, `THETA_INIT`) replace the magic numbers `9.81`, `285`, `1.5` that appeared repeatedly inline
- Added header comments explaining what the two schemes represent physically

Nothing about the numerics, grid, or time-stepping was altered, so `RunSEBBenchmark` will produce identical output to your original file.
