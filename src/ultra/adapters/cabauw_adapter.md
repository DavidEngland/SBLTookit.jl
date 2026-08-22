`CabauwAdapter` standardizes vertical tower observations from the 200m CESAR meteorological mast (Cabauw, Netherlands) into unified `MeteorologicalProfile` structures. This enables site-agnostic ingestion across both GABLS benchmark datasets and Finnish SMEAR observational networks.

**Observational Tower Comparison**

| Mast Site | Max Height ($z_{\text{max}}$) | Standard Nodes ($z$) | Displacement ($d_0$) | Roughness ($z_{0m}$) | Primary Surface |
| --- | --- | --- | --- | --- | --- |
| **Cabauw** | 200 m | 2, 10, 20, 40, 80, 140, 200 m | 0.0 m | 0.15 m | Flat Grassland |
| **Hyytiälä** | 125 m | 4.2, 8.4, 16.8, 33.6, 50.4, 67.2, 125.0 m | 12.0 m | 0.80 m | Boreal Pine Forest |
| **Värriö** | 15 m | 2, 4, 8, 15 m | 0.0 m | 0.05 m | Subarctic Fell / Forest |

---

**Cross-Site Unified Processing Pipeline**

This script ingests Cabauw DataFrame records, converts them to `StandardizedBLObservation` instances, projects the vertical temperature profile onto logarithmic Chebyshev modes (`SpectralEngine`), and evaluates stability (`UltraStability`):

```julia
using DataFrames, Dates
using .CoreTypes
using .CabauwAdapter
using .SpectralEngine
using .UltraStability

# 1. Create mock GABLS3 benchmark DataFrame row
df_cabauw = DataFrame(
    datetime = [DateTime(2026, 8, 22, 12, 0)],
    TA_2m = [18.5], TA_10m = [18.2], TA_20m = [17.8],
    TA_40m = [17.3], TA_80m = [16.5], TA_140m = [15.6], TA_200m = [14.8],
    ustar = [0.32], L_obukhov = [-50.0]
)

# 2. Extract standardized profile and observational metrics
profiles = extract_temperature_profiles(df_cabauw)
observations = extract_temperature_observations(df_cabauw; campaign="GABLS3")

prof = profiles[1]
obs = observations[1]

# 3. Spectral decomposition into 4 Chebyshev modes
cheb_coeffs = chebyshev_fingerprint(prof; n_coeffs=4, height_mapping=:log)

# 4. Compute stability regime at reference height (z_ref = 20m)
zeta = prof.metadata.reference_height / prof.metadata.obukhov_length
regime = classify_stability(zeta)

```

Would you like to build an automated site-routing interface that selects between `CabauwAdapter` and `SmearAdapter` based on input DataFrame column names?