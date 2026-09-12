Option Explicit

Public Type SEBModelConfig
    Nz As Long
    z_max As Double
    dt As Double
    t_end As Double
    Cs As Double
    Ks As Double
    T_deep As Double
    Rn0 As Double
    rho_cp As Double
    Ri_c As Double
    alpha As Double
    epsilon_c As Double
    eps_hyper As Double
    l_m As Double
    S_eff As Double
End Type

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

Private Function Phi_eps(ByVal x As Double, ByVal eps_h As Double) As Double
    Phi_eps = Sqr(x * x + eps_h * eps_h)
End Function

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

Public Sub RunSEBSimulation(ByVal use_cz0hr As Boolean, ByRef cfg As SEBModelConfig, _
                            ByRef final_Ts As Double, ByRef delta_Ts As Double, _
                            ByRef final_H0 As Double, ByRef final_Kh As Double)
    Dim Nz As Long: Nz = cfg.Nz
    Dim z() As Double: ReDim z(1 To Nz)
    Dim dz As Double: dz = (cfg.z_max - 1.5) / (Nz - 1)
    Dim g_grav As Double: g_grav = 9.81
    Dim theta_0 As Double: theta_0 = 285#

    Dim i As Long
    For i = 1 To Nz
        z(i) = 1.5 + (i - 1) * dz
    Next i

    Dim theta() As Double: ReDim theta(1 To Nz)
    Dim Ts As Double: Ts = 285#
    For i = 1 To Nz
        theta(i) = 285#
    Next i

    Dim Nt As Long: Nt = CLng(cfg.t_end / cfg.dt)
    Dim n As Long

    ReDim sub_diag(1 To Nz) As Double
    ReDim main_diag(1 To Nz) As Double
    ReDim sup_diag(1 To Nz) As Double
    ReDim rhs(1 To Nz) As Double

    Dim dtheta_dz() As Double: ReDim dtheta_dz(1 To Nz)
    Dim Ri_raw() As Double: ReDim Ri_raw(1 To Nz)
    Dim Kh() As Double: ReDim Kh(1 To Nz)
    Dim d2Ri_dz2() As Double: ReDim d2Ri_dz2(1 To Nz)
    Dim H0 As Double

    For n = 1 To Nt
        ' 1. Compute Raw Richardson Profile
        dtheta_dz(1) = (theta(1) - Ts) / z(1)
        For i = 2 To Nz
            dtheta_dz(i) = (theta(i) - theta(i - 1)) / dz
        Next i

        For i = 1 To Nz
            Ri_raw(i) = (g_grav / theta_0) * dtheta_dz(i) / (cfg.S_eff ^ 2)
        Next i

        ' 2. Evaluate Diffusivity Kh
        If Not use_cz0hr Then
            For i = 1 To Nz
                Dim g_val As Double: g_val = 1# - Ri_raw(i) / cfg.Ri_c
                If g_val < 0# Then g_val = 0#
                Kh(i) = (cfg.l_m ^ 2) * cfg.S_eff * (g_val ^ 2)
            Next i
        Else
            ' C-Z0HR C^inf Closure
            d2Ri_dz2(1) = 0#
            For i = 2 To Nz - 1
                d2Ri_dz2(i) = (Ri_raw(i + 1) - 2# * Ri_raw(i) + Ri_raw(i - 1)) / (dz ^ 2)
            Next i
            d2Ri_dz2(Nz) = d2Ri_dz2(Nz - 1)

            For i = 1 To Nz
                Dim x_val As Double: x_val = Ri_raw(i) - cfg.Ri_c
                Dim Phi As Double: Phi = Phi_eps(x_val, cfg.eps_hyper)
                Dim C_scale As Double: C_scale = 0.5 * Abs(d2Ri_dz2(i)) * (dz ^ 2) + cfg.epsilon_c * cfg.Ri_c

                Dim Ri_reg As Double
                Ri_reg = cfg.Ri_c + x_val * ((Phi + C_scale) / (Phi + (1# + cfg.alpha) * C_scale))

                Dim g_reg As Double: g_reg = 1# - Ri_reg / cfg.Ri_c
                Dim softplus_val As Double: softplus_val = 0.5 * (g_reg + Phi_eps(g_reg, cfg.eps_hyper))
                Kh(i) = (cfg.l_m ^ 2) * cfg.S_eff * (softplus_val ^ 2)
            Next i
        End If

        ' 3. SEB Surface Temperature Update
        Dim gamma_coupling As Double: gamma_coupling = cfg.rho_cp * Kh(1) / z(1)
        Dim beta_denom As Double: beta_denom = (cfg.Cs / cfg.dt) + gamma_coupling + cfg.Ks

        Ts = ((cfg.Cs / cfg.dt) * Ts + cfg.Rn0 + gamma_coupling * theta(1) + cfg.Ks * cfg.T_deep) / beta_denom

        ' Sensible heat flux: positive upward [W/m^2] (H0 < 0 indicates downward heat exchange)
        H0 = cfg.rho_cp * Kh(1) * (Ts - theta(1)) / z(1)

        ' 4. Implicit Column Diffusion (Nodal Form)
        Dim r_s As Double: r_s = Kh(1) * cfg.dt / (z(1) ^ 2)
        Dim r_v As Double: r_v = Kh(1) * cfg.dt / (dz ^ 2)

        main_diag(1) = 1# + r_s + r_v
        sup_diag(1) = -r_v
        sub_diag(1) = 0#
        rhs(1) = theta(1) + r_s * Ts

        For i = 2 To Nz - 1
            Dim r_i As Double: r_i = Kh(i) * cfg.dt / (dz ^ 2)
            sub_diag(i) = -r_i
            main_diag(i) = 1# + 2# * r_i
            sup_diag(i) = -r_i
            rhs(i) = theta(i)
        Next i

        Dim r_N As Double: r_N = Kh(Nz) * cfg.dt / (dz ^ 2)
        sub_diag(Nz) = -r_N
        main_diag(Nz) = 1# + r_N
        sup_diag(Nz) = 0#
        rhs(Nz) = theta(Nz)

        SolveTridiagonal sub_diag, main_diag, sup_diag, rhs, theta, Nz
    Next n

    final_Ts = Ts
    delta_Ts = Ts - 285#
    final_H0 = H0
    final_Kh = Kh(1)
End Sub

Public Sub RunSEBBenchmark()
    Dim grid_sizes Variant
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
        Debug.Print Format("Unregularized", "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(Ts_u, "0.00"), "@@@@@@@@@@@@") & " | " & _
                    Format(Format(dTs_u, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(H0_u, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(Kh_u, "0.00E+00"), "@@@@@@@@@@@@@@@")

        Debug.Print Format("C-Z0HR Reg.", "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(Ts_c, "0.00"), "@@@@@@@@@@@@") & " | " & _
                    Format(Format(dTs_c, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(H0_c, "0.00"), "@@@@@@@@@@@@@@@") & " | " & _
                    Format(Format(Kh_c, "0.00E+00"), "@@@@@@@@@@@@@@@")
        Debug.Print "=========================================================================" & vbCrLf
    Next idx
End Sub