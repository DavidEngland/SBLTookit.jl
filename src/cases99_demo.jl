# ==============================================================================
# Demonstration of Diagnostic Curvature Auditing
# ==============================================================================
function run_gspt_audit_demo()
    # Typical CASES-99 tower geometry (N = 7)
    z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]

    # Model parameters
    β_m, β_h = 1.0, 3.0
    L0 = 20.0

    # Nocturnal LLJ Flux Divergence Profile: L collapses near Jet axis (z ≈ 35 m)
    L_profile = [L0 * (1.0 - 0.45 * sin(π * zi / 60.0)) for zi in z_tower]

    # Exact Richardson Profile under height-varying L
    ζ_true = z_tower ./ L_profile
    Ri_obs = [Ri_model(ζ_true[i], β_m, β_h) for i in 1:length(z_tower)]

    # Run the curvature audit
    C_const, C_coord, Ri_exact, Ri_zz_obs, E_error =
        audit_gspt_curvature(z_tower, Ri_obs, L_profile, β_m, β_h)

    println("="*105)
    println("                                  GSPT CURVATURE AUDIT PROFILE")
    println("="*105)
    @printf("%-6s | %-16s | %-16s | %-16s | %-16s | %-12s\n",
            "z (m)", "Observed Curvature", "Intrinsic C_const", "Coordinate C_coord", "Exact Curvature", "Audit Error")
    println("-"*105)
    for i in 1:length(z_tower)
        @printf("%6.1f | %16.6f | %16.6f | %16.6f | %16.6f | %12.6e\n",
                z_tower[i], Ri_zz_obs[i], C_const[i], C_coord[i], Ri_exact[i], E_error[i])
    end
    println("="*105)
end

run_gspt_audit_demo()