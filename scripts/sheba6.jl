# ==============================================================================
# STAGE 6: CLOSURE COMPARISON & SCM ERROR PARTITIONING MODULE
# ==============================================================================

struct Stage6ComparisonResult
    z::Vector{Float64}
    Ri_g::Vector{Float64}
    zeta_BD::Vector{Float64}
    zeta_G07::Vector{Float64}
    delta_zeta::Vector{Float64}          # ζ_G07 - ζ_BD
    Delta_closure::Vector{Float64}       # Structural closure residual (C_const_G07 - C_const_BD)
    Delta_geom::Vector{Float64}          # Coordinate distortion residual (C_map_G07 - C_map_BD)
    tau_discretization::Vector{Float64}  # Spatial Taylor truncation correction (τ_Δz)
    fold_ratio_BD::Vector{Float64}       # |C_mapping| / (|C_const| + |C_mapping|)
    fold_ratio_G07::Vector{Float64}
end

"""
    run_stage6_closure_comparison(z, theta_z, U_z, L_profile, L_p_profile, L_pp_profile; theta_ref=273.15, g=9.81)

Executes Stage 6 comparison between Businger-Dyer and Grachev (2007) closures. Enforces
primitive-first regularization Q[M_δ u, M_δ v, M_δ θ_v] to avoid non-commutation operator artifacts.
"""
function run_stage6_closure_comparison(
    z::Vector{Float64},
    theta_z::Vector{Float64},
    U_z::Vector{Float64},
    L_profile::Vector{Float64},
    L_p_profile::Vector{Float64},
    L_pp_profile::Vector{Float64};
    theta_ref=273.15,
    g=9.81
)
    N = length(z)

    # Primitive-First Derivative Quotient: Q[M_δ θ_v, M_δ u, M_δ v]
    U_z_guarded = max.(abs.(U_z), 1e-4)
    Ri_g = (g ./ theta_ref) .* theta_z ./ (U_z_guarded .^ 2)
    Ri_g_guarded = min.(Ri_g, 0.19)

    zeta_BD = zeros(N)
    zeta_G07 = zeros(N)
    C_const_BD = zeros(N)
    C_mapping_BD = zeros(N)
    C_const_G07 = zeros(N)
    C_mapping_G07 = zeros(N)

    for i in 1:N
        # --- Businger-Dyer Closure ---
        z_bd = Ri_g_guarded[i] / (1.0 - 5.0 * Ri_g_guarded[i])
        zeta_BD[i] = z_bd
        denom_bd = 1.0 + 5.0 * z_bd
        R_zeta_bd = 1.0 / (denom_bd^2)
        R_zeta_zeta_bd = -10.0 / (denom_bd^3)

        # --- Grachev et al. (2007) Closure ---
        z_g07, dzeta_dRi = invert_grachev_Ri(Ri_g_guarded[i])
        zeta_G07[i] = z_g07
        phi_m, phi_h = eval_grachev_phi(z_g07)
        R_zeta_g07 = 1.0 / max(dzeta_dRi, 1e-6)
        R_zeta_zeta_g07 = -2.0 * 5.0 / ((1.0 + 5.0 * z_g07)^3)

        # --- Coordinate Kinematics & Mapping Terms ---
        L = L_profile[i]
        L_p = L_p_profile[i]
        L_pp = L_pp_profile[i]

        # BD Kinematics
        zeta_z_bd = (1.0 - z_bd * L_p) / L
        zeta_zz_bd = -2.0 * (L_p / L) * zeta_z_bd - (z_bd * L_pp / L)
        C_const_BD[i] = R_zeta_zeta_bd * (zeta_z_bd^2)
        C_mapping_BD[i] = R_zeta_bd * zeta_zz_bd

        # G07 Kinematics
        zeta_z_g07 = (1.0 - z_g07 * L_p) / L
        zeta_zz_g07 = -2.0 * (L_p / L) * zeta_z_g07 - (z_g07 * L_pp / L)
        C_const_G07[i] = R_zeta_zeta_g07 * (zeta_z_g07^2)
        C_mapping_G07[i] = R_zeta_g07 * zeta_zz_g07
    end

    delta_zeta = zeta_G07 .- zeta_BD
    Delta_closure = C_const_G07 .- C_const_BD
    Delta_geom = C_mapping_G07 .- C_mapping_BD

    # Discrete Taylor Truncation Correction (τ_Δz ≈ (Δz)^2 / 12 * Δ^4 Ri / Δz^4)
    dz = N > 1 ? (z[end] - z[1]) / (N - 1) : 1.0
    Ri_4th_diff = zeros(N)
    if N >= 5
        for i in 3:(N-2)
            Ri_4th_diff[i] = (Ri_g[i+2] - 4*Ri_g[i+1] + 6*Ri_g[i] - 4*Ri_g[i-1] + Ri_g[i-2]) / (dz^4)
        end
    end
    tau_discretization = ((dz^2) / 12.0) .* Ri_4th_diff

    # Fold Ratio Evaluations
    fold_ratio_BD = abs.(C_mapping_BD) ./ (abs.(C_const_BD) .+ abs.(C_mapping_BD) .+ 1e-12)
    fold_ratio_G07 = abs.(C_mapping_G07) ./ (abs.(C_const_G07) .+ abs.(C_mapping_G07) .+ 1e-12)

    return Stage6ComparisonResult(
        z, Ri_g, zeta_BD, zeta_G07, delta_zeta,
        Delta_closure, Delta_geom, tau_discretization,
        fold_ratio_BD, fold_ratio_G07
    )
end