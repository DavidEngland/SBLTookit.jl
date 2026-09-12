import numpy as np
import matplotlib.pyplot as plt

# Configuration
Ri_c = 0.20
alpha = 2.0
epsilon_c = 0.05
eps_hyper = 1e-3
eta = 1e-4

Nz = 300
z = np.linspace(1.0, 150.0, Nz)
dz = z[1] - z[0]

# Low-Level Jet Kinematics
z_J = 50.0   # Jet core height [m]
U_J = 10.0   # Jet core speed [m/s]

U = U_J * (z / z_J) * np.exp(1.0 - z / z_J)
dU_dz = (U_J / z_J) * (1.0 - z / z_J) * np.exp(1.0 - z / z_J)
S_eff = np.sqrt(dU_dz**2 + 1e-6)

g = 9.81
theta_0 = 290.0
dtheta_dz = 0.015 + 0.01 * np.exp(-z / 30.0)
N2 = (g / theta_0) * dtheta_dz

Ri_raw = N2 / (S_eff**2)

d2Ri_dz2 = np.zeros(Nz)
d2Ri_dz2[1:-1] = (Ri_raw[2:] - 2.0 * Ri_raw[1:-1] + Ri_raw[:-2]) / (dz**2)
d2Ri_dz2[0] = d2Ri_dz2[1]
d2Ri_dz2[-1] = d2Ri_dz2[-2]

# C-Z0HR mapping
x = Ri_raw - Ri_c
Phi = np.sqrt(x**2 + eps_hyper**2)
C_scale = 0.5 * np.sqrt(d2Ri_dz2**2 + eta**2) * (dz**2) + epsilon_c * Ri_c

num_reg = Phi + C_scale
den_reg = Phi + (1.0 + alpha) * C_scale
Ri_reg = Ri_c + x * (num_reg / den_reg)

l_m = 15.0 * (z / (z + 20.0))

# Exchange coefficients
g_raw = np.maximum(0.0, 1.0 - Ri_raw / Ri_c)
Sm_raw_C0 = g_raw**2
Km_raw = (l_m**2) * S_eff * Sm_raw_C0

g_reg = 1.0 - (Ri_reg / Ri_c)
# Regularized ramp R_eps
R_eps = 0.5 * (g_reg + np.sqrt(g_reg**2 + eps_hyper**2))
Sm_reg_Cinf = R_eps**2
Km_reg_Cinf = (l_m**2) * S_eff * Sm_reg_Cinf

# Create plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.8), dpi=300)

# Panel (a): Richardson number profiles
ax1.semilogx(Ri_raw, z, 'k--', label=r'Raw Richardson ($Ri_g$)', linewidth=1.5)
ax1.semilogx(Ri_reg, z, 'b-', label=r'C-Z0HR Regularized ($\widetilde{Ri}_g$)', linewidth=2.0)
ax1.axvline(Ri_c, color='r', linestyle=':', label=r'Critical Threshold ($Ri_c = 0.20$)', linewidth=1.2)
ax1.set_xlabel(r'Gradient Richardson Number, $Ri$ (log scale)', fontsize=11)
ax1.set_ylabel(r'Height, $z$ (m)', fontsize=11)
ax1.set_title(r'(a) Richardson Number Compression', fontsize=12, fontweight='bold')
ax1.grid(True, which="both", ls="--", alpha=0.5)
ax1.set_ylim(0, 150)
ax1.legend(loc='upper right', fontsize=9, framealpha=0.9)

# Panel (b): Eddy diffusivity profiles
ax2.plot(Km_raw, z, 'k--', label=r'Raw Unregularized ($K_{m,\mathrm{raw}}$)', linewidth=1.5)
ax2.plot(Km_reg_Cinf, z, 'r-', label=r'$C^\infty$ Regularized Ramp ($K_{m,\mathrm{reg}}$)', linewidth=2.0)
ax2.set_xlabel(r'Eddy Diffusivity, $K_m$ ($\mathrm{m}^2\,\mathrm{s}^{-1}$)', fontsize=11)
ax2.set_ylabel(r'Height, $z$ (m)', fontsize=11)
ax2.set_title(r'(b) Eddy Diffusivity & Residual Transport', fontsize=12, fontweight='bold')
ax2.grid(True, ls="--", alpha=0.5)
ax2.set_ylim(0, 150)
ax2.set_xlim(-0.1, max(Km_raw) * 1.05)
ax2.legend(loc='upper right', fontsize=9, framealpha=0.9)

# Inset for panel (b) showing log-scale supercritical tail
ax2_inset = ax2.inset_axes([0.45, 0.25, 0.50, 0.40])
ax2_inset.semilogx(Km_reg_Cinf[z > 50], z[z > 50], 'r-', linewidth=1.5)
ax2_inset.set_xlabel(r'$K_m$ (log scale)', fontsize=8)
ax2_inset.set_ylabel(r'$z$ (m)', fontsize=8)
ax2_inset.tick_params(labelsize=7)
ax2_inset.set_title(r'$C^\infty$ Supercritical Tail', fontsize=8, fontweight='bold')
ax2_inset.grid(True, which="both", ls=":", alpha=0.6)

plt.tight_layout()
plt.savefig('cz0hr_llj_diagnostics.png', dpi=300)
print("Figure generated successfully.")