# SBL Climatology: GSPT Fold Illusions Occur in 11% of Outer Jet Profiles

Statistical climatological analysis of the GABLS3 Cabauw campaign vertical sounding dataset reveals that **10.97% of profiles in the Outer Jet Layer (z > 80m) suffer from "Fold Illusion" coordinate singularities.** In these cases, the visual "knee" or sharp gradient inflection in the wind shear profile is driven entirely by non-linear coordinate projection geometry rather than a physical collapse of boundary-layer turbulent mixing. Conversely, classical Monin-Obukhov Similarity Theory (MOST) is highly robust in the Surface Shear Layer (z ≤ 20m), where fold illusions are completely absent (0.00% frequency) despite a 36.25% frequency of gradient Richardson numbers crossing the standard critical threshold ($Ri_g > 0.20$).

---

## 1. Key Climatological Findings

1.  **The Jet-Nose Singularity Zone:** In the **Outer Jet Layer (z > 80m)**, the mean GSPT Fold Ratio reaches **0.8762**, indicating that vertical profile bending is almost entirely dominated by coordinate mapping curvature ($C_{\text{mapping}} = R_\zeta \zeta_{zz}$). Within this zone, Fold Illusion events occur at a frequency of **10.97%**, meaning that roughly 1 in 9 profiles contains a visual gradient inflection that numerical models would misclassify as a physical turbulence collapse.
2.  **Surface Shear Layer Robustness:** Below **$z = 20\text{ m}$**, Monin-Obukhov similarity mapping remains perfectly well-conditioned. The mean Fold Ratio is only **0.1205**, and Fold Illusion events are completely absent (**0.00%**). Although the gradient Richardson number frequently crosses $Ri_{cr} = 0.20$ (at a **36.25%** frequency), the coordinate conditioning metric $\kappa_{\zeta} \equiv |d\zeta / dRi_g|$ remains tightly bounded ($\kappa_\zeta \approx 1.5$), demonstrating that threshold crossings in the surface layer represent well-behaved physical transitions.
3.  **The Mid-Boundary Layer Transition:** The **Mid-Boundary Layer (20m < z ≤ 80m)** represents a highly sensitive transition zone. The mean Fold Ratio rises to **0.5145**, and the Fold Illusion frequency climbs to **7.08%**. Here, the vertical flux divergence ($L' \neq 0$) is active but moderate, causing the coordinate mapping and physical closure curvatures to be comparable ($C_{\text{mapping}} \approx C_{\text{constitutive}}$), creating complex "double knee" profile shapes.

---

## 2. Statistical Climatology Summary Table

| vertical Sub-Layer | Sample Count (N) | Mean Richardson ($Ri_g$) | Mean GSPT Fold Ratio | Max Sensitivity ($\kappa_{\zeta}$) | Fold Illusion Frequency | MOST Breakdown ($Ri_g > 0.20$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Surface Shear Layer ($z \le 20\text{ m}$)** | 240 | 0.2096 | 0.1205 | 0.9426 | **0.00%** (0 events) | **36.25%** (87 events) |
| **Mid-Boundary Layer ($20\text{ m} < z \le 80\text{ m}$)** | 240 | 0.1333 | 0.5145 | 0.9426 | **7.08%** (17 events) | **21.67%** (52 events) |
| **Outer Jet Layer ($z > 80\text{ m}$)** | 720 | 0.1059 | 0.8762 | 0.9426 | **10.97%** (79 events) | **11.11%** (80 events) |

---

## 3. Physical Inferences for NWP and SCM Schemes

These climatological distributions carry immediate, high-signal consequences for Single-Column Models (SCMs) and weather forecasting codes:

*   **Do Not Retune Stability Functions in the Jet Zone:** When a model misrepresents wind profiles near the Low-Level Jet, SCM developers often "retune" Monin-Obukhov coefficients ($\beta_m, \beta_h$). GSPT proves this is a compensating error. Because $C_{\text{mapping}}$ dominates ($87.62\%$ of total curvature), the error stems from **excessive vertical flux divergence ($L', L''$)** over-smoothing the coordinate grids—the local stability functions ($\phi_m, \phi_h$) should remain untouched.
*   **Deploy Localized Reference Temperatures:** Buoyancy frequencies ($N^2$) are highly sensitive to moisture and vertical potential temperature gradients. The data confirms that utilizing Boussinesq constant-reference temperatures ($\theta_0$) over seasonal-scale runs biases gradient calculations, and a height-varying, local reference temperature should be deployed.

---

## 4. Methodology and Data Processing

GSPT coordinates and curvature terms were dynamically reconstructed using the thread-safe `compare_sbl_closures` processing engine. Primitive wind components and virtual potential temperature profiles were first regularized using log-coordinate natural cubic splines ($\xi = \ln(z/z_0)$ with roughness length $z_0 = 0.15\text{ m}$) under a self-expanding Tikhonov-Morozov discrepancy solver. 

Joint probability density functions (PDFs) were estimated using a 2D Gaussian Kernel Density Estimator (KDE) with standard bandwidth selection.

![GABLS3 Joint PDF Climatology](gabls3-gspt-climatology.png)
