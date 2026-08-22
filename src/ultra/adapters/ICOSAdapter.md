`ICOSAdapter.jl` bridges sparse two-level tower installations (such as standard ICOS greenhouse gas stations) to multi-level spectral engines by using Monin-Obukhov Similarity Theory (MOST) to reconstruct physically consistent boundary layer profiles.

**MOST Stability Coordinate Formulations ($\chi$)**

| Regime | Condition ($\zeta = z/L$) | Integrated Correction $\psi(\zeta)$ | Modified Coordinate $\chi(z)$ |
| --- | --- | --- | --- |
| **Neutral** | $\zeta \to 0$ or `:none` | $0$ | $\ln z$ |
| **Stable** ($\psi_m, \psi_h$) | $\zeta > 0$ | $-5\zeta$ | $\ln z + 5(z/L)$ |
| **Unstable Momentum** ($\psi_m$) | $\zeta \le 0$ | $2\ln\left(\frac{1+x}{2}\right) + \ln\left(\frac{1+x^2}{2}\right) - 2\arctan(x) + \frac{\pi}{2}$ | $\ln z - \psi_m(z/L)$ |
| **Unstable Heat** ($\psi_h$) | $\zeta \le 0$ | $2\ln\left(\frac{1+y}{2}\right)$ | $\ln z - \psi_h(z/L)$ |

*Where $x = (1 - 16\zeta)^{1/4}$ and $y = (1 - 16\zeta)^{1/2}$.*

---

**End-to-End Pipeline: Sparse ICOS Upscaling to Chebyshev Decomposition**

This pipeline upscales a two-level observation (10m and 50m) into a five-level virtual profile to satisfy the $N \ge 3$ level constraint required for `SpectralEngine` modal decomposition:

```julia
using Dates
using .CoreTypes
using .ICOSAdapter
using .SpectralEngine
using .UltraStability

# 1. Sparse 2-level ICOS temperature measurement
dt = DateTime(2026, 8, 22, 14, 0)
z_low, T_low   = 10.0, 18.5  # 10m temperature (°C)
z_high, T_high = 50.0, 16.2  # 50m temperature (°C)
ustar, L_obukhov = 0.25, -35.0 # Convective conditions (L < 0)

# 2. Upscale 2-level observation to 5 virtual target levels using MOST heat scaling (:psi_h)
obs = upscale_sparse_icos_observation(
    dt, z_low, z_high, T_low, T_high, ustar, L_obukhov;
    campaign="ICOS-DE-S3",
    target_levels=[2.0, 5.0, 10.0, 20.0, 50.0],
    stability_correction=:psi_h
)

# 3. Construct MeteorologicalProfile and extract 4-mode Chebyshev fingerprint
meta = ProfileMetadata(dt, ustar, L_obukhov, z_high, 0.0)
prof = MeteorologicalProfile(meta, obs.heights, obs.values)
c_modes = chebyshev_fingerprint(prof; n_coeffs=4, height_mapping=:log)

# 4. Classify stability state
regime = classify_stability(z_low / L_obukhov)

```

Would you like to write a unified multi-adapter pipeline (`UnifiedBLIngestion.jl`) that auto-detects input file formats (SMEAR, Cabauw, AmeriFlux, or ICOS) and standardizes them into `StandardizedBLObservation` objects?