Option Explicit

'=============================================================================
' C-Z0HR REGULARIZATION BENCHMARK vs. UNREGULARIZED CLOSURE IN A LOW-LEVEL JET
' SBLToolkit.jl / Generalized Similarity Profile Theory (GSPT)
' Macro and User-Defined Function (UDF) Suite for Microsoft Excel / VBA
'=============================================================================
' Core Theoretical Principles:
' 1. C^∞ Hyperbolic Metric:   Φ_ε(x) = √(x² + ε_hyper²)
' 2. Curvature Scale:          C(z, Δz) = 0.5 * |∂²Ri_g/∂z²| * (Δz)² + ε_c * Ri_c
' 3. C-Z0HR Regularized Ri:    Ri_reg = Ri_c + x * (Φ + C) / (Φ + (1 + α) * C)
' 4. Partial Sensitivity:      D_α^(R) ≡ ∂Ri_g^reg / ∂Ri_g |_C
'    Threshold Crossing:       D_α^(R)(Ri_c) ≈ 1 / (1 + α) = 0.333 (for α = 2.0)
' 5. Smooth C^∞ Diffusivity:   Sm_reg = [0.5 * (g + √(g² + ε_hyper²))]²
'                              where g = 1 - Ri_reg / Ri_c
'=============================================================================

'-----------------------------------------------------------------------------
' PUBLIC USER-DEFINED FUNCTIONS (UDFs) FOR WORKSHEET FORMULAS
'-----------------------------------------------------------------------------

''' <summary>
''' Smooth C^∞ hyperbolic distance metric replacing non-differentiable |x| kinks.
''' </summary>
Public Function CZ0HR_Phi(ByVal x As Double, Optional ByVal eps_hyper As Double = 0.001) As Double
    CZ0HR_Phi = Sqr(x * x + eps_hyper * eps_hyper)
End Function

''' <summary>
''' Evaluates Curvature-Aware Zero-Offset Hyperbolic Regularization (C-Z0HR) on Ri_g.
''' </summary>
Public Function CZ0HR_Regulate( _
    ByVal Ri_raw As Double, _
    ByVal d2Ri_dz2 As Double, _
    ByVal dz As Double, _
    Optional ByVal Ri_c As Double = 0.2, _
    Optional ByVal alpha_param As Double = 2#, _
    Optional ByVal epsilon_c As Double = 0.05, _
    Optional ByVal eps_hyper As Double = 0.001 _
) As Double
    Dim x As Double, Phi As Double, C As Double
    Dim num_reg As Double, den_reg As Double

    x = Ri_raw - Ri_c
    Phi = Sqr(x * x + eps_hyper * eps_hyper)
    C = 0.5 * Abs(d2Ri_dz2) * (dz * dz) + (epsilon_c * Ri_c)

    num_reg = Phi + C
    den_reg = Phi + (1# + alpha_param) * C

    CZ0HR_Regulate = Ri_c + x * (num_reg / den_reg)
End Function

''' <summary>
''' Computes the partial sensitivity derivative D_α^(R) = ∂Ri_reg / ∂Ri_raw |_C.
''' </summary>
Public Function CZ0HR_DAlpha( _
    ByVal Ri_raw As Double, _
    ByVal d2Ri_dz2 As Double, _
    ByVal dz As Double, _
    Optional ByVal Ri_c As Double = 0.2, _
    Optional ByVal alpha_param As Double = 2#, _
    Optional ByVal epsilon_c As Double = 0.05, _
    Optional ByVal eps_hyper As Double = 0.001 _
) As Double
    Dim x As Double, Phi As Double, C As Double
    Dim num_D As Double, den_D As Double

    x = Ri_raw - Ri_c
    Phi = Sqr(x * x + eps_hyper * eps_hyper)
    C = 0.5 * Abs(d2Ri_dz2) * (dz * dz) + (epsilon_c * Ri_c)

    num_D = (Phi * Phi) + 2# * (1# + alpha_param) * C * Phi + (1# + alpha_param) * (C * C)
    den_D = (Phi + (1# + alpha_param) * C) ^ 2

    CZ0HR_DAlpha = num_D / den_D
End Function

''' <summary>
''' Classical piecewise C^0 short-tail stability function Sm(Ri) = max(0, 1 - Ri/Ri_c)².
''' </summary>
Public Function Sm_ShortTail_C0(ByVal Ri As Double, Optional ByVal Ri_c As Double = 0.2) As Double
    Dim val As Double
    val = 1# - (Ri / Ri_c)
    If val > 0# Then
        Sm_ShortTail_C0 = val * val
    Else
        Sm_ShortTail_C0 = 0#
    End If
End Function

''' <summary>
''' C^∞ smoothed short-tail stability function using hyperbolic activation.
''' </summary>
Public Function Sm_ShortTail_Cinf( _
    ByVal Ri As Double, _
    Optional ByVal Ri_c As Double = 0.2, _
    Optional ByVal eps_hyper As Double = 0.001 _
) As Double
    Dim g As Double, g_smooth As Double
    g = 1# - (Ri / Ri_c)
    g_smooth = 0.5 * (g + Sqr(g * g + eps_hyper * eps_hyper))
    Sm_ShortTail_Cinf = g_smooth * g_smooth
End Function

'-----------------------------------------------------------------------------
' BENCHMARK EXECUTION SUBROUTINE
'-----------------------------------------------------------------------------

''' <summary>
''' Simulates a 1D Low-Level Jet (LLJ) vertical profile (0–150 m) and audits
''' C-Z0HR singularity suppression and eddy diffusivity recovery.
''' Writes complete profile data and summary metrics to the ActiveSheet.
''' </summary>
Public Sub Run_CZ0HR_Benchmark()
    On Error GoTo ErrorHandler

    ' 1. Grid & Physical Constants
    Const Nz As Long = 300
    Const z_min As Double = 1#
    Const z_max As Double = 150#
    Const z_J As Double = 50#          ' Jet core height [m]
    Const U_J As Double = 10#          ' Jet core wind speed [m/s]
    Const g As Double = 9.81           ' Gravitational acceleration [m/s²]
    Const theta_0 As Double = 290#     ' Reference potential temperature [K]

    ' 2. Regularization Hyperparameters
    Const Ri_c As Double = 0.2         ' Critical Richardson threshold
    Const alpha_param As Double = 2#   ' Hyperbolic damping control parameter
    Const epsilon_c As Double = 0.05   ' Curvature scale cap factor
    Const eps_hyper As Double = 0.001  ' C^∞ hyperbolic smoothing parameter

    Dim dz As Double
    dz = (z_max - z_min) / (Nz - 1)

    ' 3. Allocate Computation Arrays
    Dim z(1 To Nz) As Double, U(1 To Nz) As Double, dU_dz(1 To Nz) As Double
    Dim S_eff(1 To Nz) As Double, dtheta_dz(1 To Nz) As Double, N2(1 To Nz) As Double
    Dim Ri_raw(1 To Nz) As Double, d2Ri_dz2(1 To Nz) As Double
    Dim Ri_reg(1 To Nz) As Double, D_alpha(1 To Nz) As Double, C_scale(1 To Nz) As Double
    Dim l_m(1 To Nz) As Double, Sm_raw_C0(1 To Nz) As Double
    Dim Sm_reg_C0(1 To Nz) As Double, Sm_reg_Cinf(1 To Nz) As Double
    Dim Km_raw(1 To Nz) As Double, Km_reg_C0(1 To Nz) As Double, Km_reg_Cinf(1 To Nz) As Double

    Dim i As Long, zi As Double
    Dim x As Double, Phi As Double, C As Double
    Dim num_reg As Double, den_reg As Double, num_D As Double, den_D As Double
    Dim g_factor As Double, g_smooth As Double

    ' 4. Step 1: Kinematics & Stratification Profiles
    For i = 1 To Nz
        zi = z_min + (i - 1) * dz
        z(i) = zi

        ' Low-Level Jet wind speed & vertical shear
        U(i) = U_J * (zi / z_J) * Exp(1# - (zi / z_J))
        dU_dz(i) = (U_J / z_J) * (1# - (zi / z_J)) * Exp(1# - (zi / z_J))
        S_eff(i) = Sqr(dU_dz(i) * dU_dz(i) + 0.000001) ' 1e-6 floor prevents /0

        ' Potential temperature gradient and Brunt-Väisälä frequency squared N²
        dtheta_dz(i) = 0.015 + 0.01 * Exp(-zi / 30#)
        N2(i) = (g / theta_0) * dtheta_dz(i)

        ' Raw gradient Richardson number
        Ri_raw(i) = N2(i) / (S_eff(i) * S_eff(i))
    Next i

    ' 5. Step 2: Second Spatial Derivative (Central Finite Difference)
    For i = 2 To Nz - 1
        d2Ri_dz2(i) = (Ri_raw(i + 1) - 2# * Ri_raw(i) + Ri_raw(i - 1)) / (dz * dz)
    Next i
    d2Ri_dz2(1) = d2Ri_dz2(2)
    d2Ri_dz2(Nz) = d2Ri_dz2(Nz - 1)

    ' 6. Step 3: C-Z0HR Regularization Mapping
    For i = 1 To Nz
        x = Ri_raw(i) - Ri_c
        Phi = Sqr(x * x + eps_hyper * eps_hyper)
        C = 0.5 * Abs(d2Ri_dz2(i)) * (dz * dz) + (epsilon_c * Ri_c)
        C_scale(i) = C

        num_reg = Phi + C
        den_reg = Phi + (1# + alpha_param) * C
        Ri_reg(i) = Ri_c + x * (num_reg / den_reg)

        num_D = (Phi * Phi) + 2# * (1# + alpha_param) * C * Phi + (1# + alpha_param) * (C * C)
        den_D = den_reg * den_reg
        D_alpha(i) = num_D / den_D
    Next i

    ' 7. Step 4: Turbulent Length Scale, Stability Functions & Eddy Diffusivity Km
    For i = 1 To Nz
        ' Blackadar-type asymptotic mixing length
        l_m(i) = 15# * (z(i) / (z(i) + 20#))

        ' Classical C^0 Short-Tail Stability Function
        Sm_raw_C0(i) = Sm_ShortTail_C0(Ri_raw(i), Ri_c)
        Sm_reg_C0(i) = Sm_ShortTail_C0(Ri_reg(i), Ri_c)

        ' C^∞ Hyperbolic Smoothed Stability Function
        g_factor = 1# - (Ri_reg(i) / Ri_c)
        g_smooth = 0.5 * (g_factor + Sqr(g_factor * g_factor + eps_hyper * eps_hyper))
        Sm_reg_Cinf(i) = g_smooth * g_smooth

        ' Eddy Diffusivity Km = l_m² * S * Sm
        Km_raw(i) = (l_m(i) * l_m(i)) * S_eff(i) * Sm_raw_C0(i)
        Km_reg_C0(i) = (l_m(i) * l_m(i)) * S_eff(i) * Sm_reg_C0(i)
        Km_reg_Cinf(i) = (l_m(i) * l_m(i)) * S_eff(i) * Sm_reg_Cinf(i)
    Next i

    ' 8. Step 5: Fast Bulk Excel Output using 2D Array
    Dim prevScreenUpdating As Boolean
    Dim prevCalculation As XlCalculation
    prevScreenUpdating = Application.ScreenUpdating
    prevCalculation = Application.Calculation

    Application.ScreenUpdating = False
    Application.Calculation = xlCalculationManual

    Dim ws As Worksheet
    Set ws = ActiveSheet

    ws.Cells.Clear

    ' Title Block
    ws.Range("A1").Value = "C-Z0HR REGULARIZATION BENCHMARK AUDIT: LOW-LEVEL JET (z_J = 50 m)"
    ws.Range("A1").Font.Size = 14
    ws.Range("A1").Font.Bold = True

    ws.Range("A2").Value = "Parameters: Ri_c = " & Ri_c & " | alpha = " & alpha_param & _
                           " | eps_c = " & epsilon_c & " | eps_hyper = " & eps_hyper & _
                           " | dz = " & Format$(dz, "0.000") & " m"
    ws.Range("A2").Font.Italic = True

    ' Table Headers
    Dim headers As Variant
    headers = Array("z (m)", "U (m/s)", "S_eff (s⁻¹)", "Ri_raw", "Ri_reg", "D_alpha^(R)", _
                    "C_scale", "Km_raw (m²/s)", "Km_reg_C0 (m²/s)", "Km_reg_Cinf (m²/s)")

    ws.Range("A4:J4").Value = headers
    ws.Range("A4:J4").Font.Bold = True
    ws.Range("A4:J4").Interior.Color = RGB(220, 230, 242)

    ' Assemble Output Matrix (Nz rows x 10 cols)
    Dim outData(1 To Nz, 1 To 10) As Variant
    For i = 1 To Nz
        outData(i, 1) = z(i)
        outData(i, 2) = U(i)
        outData(i, 3) = S_eff(i)
        outData(i, 4) = Ri_raw(i)
        outData(i, 5) = Ri_reg(i)
        outData(i, 6) = D_alpha(i)
        outData(i, 7) = C_scale(i)
        outData(i, 8) = Km_raw(i)
        outData(i, 9) = Km_reg_C0(i)
        outData(i, 10) = Km_reg_Cinf(i)
    Next i

    ws.Range("A5").Resize(Nz, 10).Value = outData

    ' Number Formatting
    ws.Range("A5").Resize(Nz, 3).NumberFormat = "0.000"
    ws.Range("D5").Resize(Nz, 1).NumberFormat = "0.0000E+00"
    ws.Range("E5").Resize(Nz, 3).NumberFormat = "0.00000"
    ws.Range("H5").Resize(Nz, 3).NumberFormat = "0.0000E+00"

    ws.Columns("A:J").AutoFit

    ' Restore Excel Environment
    Application.Calculation = prevCalculation
    Application.ScreenUpdating = prevScreenUpdating

    ' 9. Immediate Window Summary
    Debug.Print "========================================================================================="
    Debug.Print "C-Z0HR BENCHMARK AUDIT COMPLETE: 300 levels generated."
    Debug.Print "Maximum Ri_raw: " & Format$(Ri_raw(100), "0.0000E+00")
    Debug.Print "Theoretical Threshold D_alpha bound: " & Format$(1# / (1# + alpha_param), "0.000")
    Debug.Print "========================================================================================="

    Exit Sub

ErrorHandler:
    Application.Calculation = xlCalculationAutomatic
    Application.ScreenUpdating = True
    MsgBox "Error during C-Z0HR benchmark: " & Err.Description, vbCritical, "C-Z0HR Error"
End Sub