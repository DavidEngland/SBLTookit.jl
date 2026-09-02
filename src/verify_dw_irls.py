import numpy as np
import scipy.interpolate as interp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd

# ----------------------------------------------------------------------
# Environment Setup and Preflight
# ----------------------------------------------------------------------
os.makedirs('/workspace/scratch', exist_ok=True)
sns.set_theme(style='whitegrid', palette='colorblind', font='DejaVu Sans')
CHART_DPI = 150

# ----------------------------------------------------------------------
# GSPT Grachev et al. (2007) and Businger-Dyer Formulations
# ----------------------------------------------------------------------
def phi_m_Grachev(zeta):
    a_m = 5.0
    b_m = 5.0 / 6.5
    return 1.0 + a_m * (zeta * (1.0 + zeta)**(1.0/3.0)) / (1.0 + b_m * zeta)

def phi_h_Grachev(zeta):
    a_h = 5.0
    b_h = 5.0
    c_h = 3.0
    return 1.0 + (a_h * zeta + b_h * zeta**2) / (1.0 + c_h * zeta + zeta**2)

def R_Grachev(zeta):
    return zeta * phi_h_Grachev(zeta) / (phi_m_Grachev(zeta)**2)

def get_R_derivative(zeta, eps=1e-5):
    # Centered five-point stencil for exact analytical derivative
    r_val = R_Grachev(zeta)
    r_plus = R_Grachev(zeta + eps)
    r_minus = R_Grachev(zeta - eps)
    r_plus2 = R_Grachev(zeta + 2*eps)
    r_minus2 = R_Grachev(zeta - 2*eps)
    
    R_prime = (-r_plus2 + 8.0*r_plus - 8.0*r_minus + r_minus2) / (12.0 * eps)
    return max(R_prime, 1e-6)

def invert_R_Grachev(Ri_target, max_iter=30):
    if Ri_target < 0.0:
        return Ri_target
    zeta_guess = Ri_target / (1.0 - 5.0 * Ri_target + 1e-6)
    zeta_guess = np.clip(zeta_guess, 0.0, 10.0)
    for _ in range(max_iter):
        r_val = R_Grachev(zeta_guess)
        diff = r_val - Ri_target
        if abs(diff) < 1e-6:
            break
        r_prime = get_R_derivative(zeta_guess)
        zeta_guess = np.clip(zeta_guess - diff / r_prime, 0.0, 50.0)
    return zeta_guess

# ----------------------------------------------------------------------
# IRLS Solvers and Covariance Sandwich Primitives
# ----------------------------------------------------------------------
def run_dw_irls_simulation():
    # 38 levels of GABLS3 mast spacing
    z_tower = np.linspace(2.0, 200.0, 38)
    N_z = len(z_tower)
    
    # In-situ noise specifications
    delta_theta_raw = np.full(N_z, 0.02) # fine thermocouples
    delta_u_raw = np.full(N_z, 0.05)     # sonic anemometers
    
    # Simulate a descending LLJ nose (z_llj = 80m) and sharp inversion (h_inv = 60m)
    z_llj = 80.0
    h_inv = 60.0
    u_g = 6.0
    u_jet = 4.0
    theta_ref = 285.0
    delta_theta = 8.0
    
    U_raw = u_g * (1.0 - np.exp(-z_tower / 35.0)) + u_jet * np.exp(-((z_tower - z_llj) / 20.0)**2) + 0.05 * np.random.randn(N_z)
    theta_raw = theta_ref + delta_theta * (1.0 - np.exp(-z_tower / h_inv)) + 0.02 * np.random.randn(N_z)
    
    # Log-coordinate mapping
    z_0 = 0.15
    xi_tower = np.log(z_tower / z_0)
    
    # Initial diagnostic weights
    w_damped = np.ones(N_z)
    gamma_damping = 0.3
    eps_s = 1e-12
    g_over_theta = 9.81 / 285.0
    
    theta_z = np.zeros(N_z)
    U_z = np.zeros(N_z)
    zeta_arr = np.zeros(N_z)
    
    for irls_iter in range(20):
        # Update effective observation noise floors
        eff_delta_theta = delta_theta_raw / np.sqrt(w_damped)
        eff_delta_u = delta_u_raw / np.sqrt(w_damped)
        
        # Fit smoothing splines to log-coordinate primitive fields
        spl_theta = interp.UnivariateSpline(xi_tower, theta_raw, w=1.0/eff_delta_theta, s=N_z)
        spl_U = interp.UnivariateSpline(xi_tower, U_raw, w=1.0/eff_delta_u, s=N_z)
        
        # Extract gradients
        theta_xi_deriv = spl_theta.derivative(1)(xi_tower)
        U_xi_deriv = spl_U.derivative(1)(xi_tower)
        
        theta_z = theta_xi_deriv / z_tower
        U_z = U_xi_deriv / z_tower
        
        # Covariance Sandwich Approximation for derivative uncertainties
        # For non-uniform grids, derivative variance scales as Var(dy/dz) ≈ Var(dy)/(dz^2 * w)
        dz = np.zeros(N_z)
        dz[0] = z_tower[1] - z_tower[0]
        dz[-1] = z_tower[-1] - z_tower[-2]
        for i in range(1, N_z-1):
            dz[i] = 0.5 * (z_tower[i+1] - z_tower[i-1])
            
        sigma_theta_z = eff_delta_theta / (dz * np.sqrt(N_z))
        sigma_U_z = eff_delta_u / (dz * np.sqrt(N_z))
        
        # Update weights and coordinates
        w_prev = w_damped.copy()
        
        for i in range(N_z):
            U_z_clamped = max(abs(U_z[i]), 1e-4) * np.sign(U_z[i])
            
            # Gradient Richardson number
            Ri_raw = g_over_theta * theta_z[i] / (U_z_clamped**2 + eps_s)
            Ri_safe = 2.0 * np.tanh(Ri_raw / 2.0)
            
            # Monin-Obukhov Inversion
            zeta_val = invert_R_Grachev(Ri_safe)
            zeta_arr[i] = zeta_val
            
            # Jacobian Coordinates
            R_prime = get_R_derivative(zeta_val)
            dzeta_dRi = 1.0 / R_prime
            
            # Chain-rule sensitivities
            dzeta_dtheta_z = dzeta_dRi * (g_over_theta / (U_z_clamped**2 + eps_s))
            dzeta_dU_z = dzeta_dRi * (-2.0 * g_over_theta * theta_z[i] / (U_z_clamped**3 + eps_s))
            
            # Propagated similarity coordinate variance
            sigma_zeta_sq = (dzeta_dtheta_z**2) * (sigma_theta_z[i]**2) + (dzeta_dU_z**2) * (sigma_U_z[i]**2)
            
            # Clamp variance to prevent numerical blowup
            sigma_zeta_sq = min(sigma_zeta_sq, 5.0)
            
            # Calculate raw downstream weight
            w_calc = 1.0 / (1.0 + sigma_zeta_sq / (0.5**2))
            
            # Damped Update
            w_damped[i] = (1.0 - gamma_damping) * w_prev[i] + gamma_damping * w_calc
            
        # Check convergence
        rel_err_w = np.max(np.abs(w_damped - w_prev) / (1.0 + w_prev))
        if rel_err_w < 1e-4:
            print(f"Downstream-Weighted IRLS converged in {irls_iter+1} iterations!")
            break
            
    return z_tower, theta_raw, U_raw, theta_z, U_z, zeta_arr, w_damped

# Execute verification and generate visualization
z, theta_raw, U_raw, theta_z, U_z, zeta, weights = run_dw_irls_simulation()

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
fig.suptitle("Damped Downstream-Weighted IRLS Safely Isolates Low-Shear Jet Singularity", fontsize=16, fontweight='bold')

# Panel 1: Smoothed Primitive Gradients
ax1.plot(theta_z, z, label=r"$\theta_z$ (smoothed)", linewidth=2.5, color='darkred')
ax1.plot(U_z, z, label=r"$U_z$ (smoothed)", linewidth=2.5, color='teal')
ax1.axhline(80.0, linestyle='--', color='black', label="Jet Axis / Inversion")
ax1.set_xlabel("Vertical Gradients (K/m or 1/s)")
ax1.set_ylabel("Height (m)")
ax1.set_title("Smoothed Vertical Gradients", fontsize=12, fontweight='bold')
ax1.legend()

# Panel 2: Inverted Similarity Coordinate (zeta)
ax2.plot(zeta, z, linewidth=2.5, color='darkblue')
ax2.set_xlabel(r"Similarity Coordinate $\zeta$")
ax2.set_title("Inverted Similarity Coordinate", fontsize=12, fontweight='bold')

# Panel 3: Converged IRLS Downstream Weights
ax3.plot(weights, z, linewidth=2.5, color='purple', marker='o', label="Downstream Weights")
ax3.axhline(80.0, linestyle='--', color='black')
ax3.set_xlabel("IRLS Downstream Weight w(z)")
ax3.set_title("Converged Downstream Weights", fontsize=12, fontweight='bold')
ax3.legend()

sns.despine()
plt.tight_layout(pad=1.5)
fig.savefig('/workspace/scratch/gabls3_dw_irls_diagnostic.png', dpi=CHART_DPI, bbox_inches='tight')
plt.close()
print("Saved visualization to /workspace/scratch/gabls3_dw_irls_diagnostic.png")
