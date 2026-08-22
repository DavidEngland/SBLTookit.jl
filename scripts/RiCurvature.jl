using Plots

# Initialize GR backend for rendering
gr()

# Common physical constants
L0 = 50.0          # Surface Monin-Obukhov length [m]
beta_m = 4.7       # Empirical stability constant

# ==============================================================================
# CASE 1: Linear Flux Divergence L(z) = L0 + L'z (Asymptotic Compression)
# ==============================================================================
z1 = range(0.1, 300.0, length=500)
L_prime1 = 0.5

L1 = L0 .+ L_prime1 .* z1
zeta1 = z1 ./ L1

# Coordinate derivatives (L'' = 0)
zeta_z1 = (1.0 .- zeta1 .* L_prime1) ./ L1
zeta_zz1 = -2.0 .* L_prime1 .* zeta_z1 ./ L1

# Similarity curvature components
Ri_zeta1 = 1.0 ./ (1.0 .+ beta_m .* zeta1) .^ 2
Ri_zetazeta1 = -2.0 .* beta_m ./ (1.0 .+ beta_m .* zeta1) .^ 3

Ri_zz_intrinsic1 = Ri_zetazeta1 .* (zeta_z1 .^ 2)
Ri_zz_coord1 = Ri_zeta1 .* zeta_zz1
Ri_zz_total1 = Ri_zz_intrinsic1 .+ Ri_zz_coord1

p1 = plot(
    [Ri_zz_intrinsic1 Ri_zz_coord1 Ri_zz_total1],
    z1,
    label=["Intrinsic (Ri_ζζ ζ_z²)" "Coordinate (Ri_ζ ζ_zz)" "Total (d²Ri/dz²)"],
    xlabel="Curvature [m⁻²]",
    ylabel="Physical Height z [m]",
    title="Case 1: Linear L(z)\n(Ratio |Coord/Int| → 10.6% at surface)",
    lw=2,
    linestyle=[:dash :dot :solid],
    legend=:topleft
)

# ==============================================================================
# CASE 2: Non-Linear Flux Divergence L(z) = L0 + c*z² (Transversality Loss)
# ==============================================================================
c = 0.005
z_star = sqrt(L0 / c) # Transversality condition z* = 100m
z2 = range(0.1, 150.0, length=500)

L2 = L0 .+ c .* (z2 .^ 2)
L_p2 = 2.0 .* c .* z2
L_pp2 = fill(2.0 * c, length(z2))

zeta2 = z2 ./ L2
zeta_z2 = (1.0 .- zeta2 .* L_p2) ./ L2
zeta_zz2 = -2.0 .* L_p2 ./ L2 .* zeta_z2 .- (zeta2 .* L_pp2 ./ L2)

# Similarity curvature components
Ri_zeta2 = 1.0 ./ (1.0 .+ beta_m .* zeta2) .^ 2
Ri_zetazeta2 = -2.0 .* beta_m ./ (1.0 .+ beta_m .* zeta2) .^ 3

Ri_zz_intrinsic2 = Ri_zetazeta2 .* (zeta_z2 .^ 2)
Ri_zz_coord2 = Ri_zeta2 .* zeta_zz2
Ri_zz_total2 = Ri_zz_intrinsic2 .+ Ri_zz_coord2

p2 = plot(
    [Ri_zz_intrinsic2 Ri_zz_coord2 Ri_zz_total2],
    z2,
    label=["Intrinsic" "Coordinate" "Total"],
    xlabel="Curvature [m⁻²]",
    ylabel="Physical Height z [m]",
    title="Case 2: Non-Linear L(z)\n(Coordinate Curvature Dominates at z* = 100m)",
    lw=2,
    linestyle=[:dash :dot :solid],
    legend=:topleft
)
hline!(p2, [z_star], color=:red, linestyle=:dashdot, label="Transversality Level z*")

# ==============================================================================
# Combine & Save Output
# ==============================================================================
plt = plot(p1, p2, layout=(1, 2), size=(1100, 500), margin=5Plots.mm)

# Display interactive plot window (if in IDE/REPL)
display(plt)

# Save visual graphic file to working directory
savefig(plt, "gspt_curvature_decomposition.png")
println("Plot figure written to 'gspt_curvature_decomposition.png'")