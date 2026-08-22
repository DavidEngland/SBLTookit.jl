The `UnifiedBLIngestion` module auto-detects observational network column signatures (Cabauw, AmeriFlux, NEON, ICOS, SMEAR) and routes them through their respective adapters into standardized boundary layer structures.

**Network Signature Detection Matrix**

| Network Format | Column Signature Criteria | Ingestion Strategy |
| --- | --- | --- |
| **`:cabauw`** | Contains `"TA_2m"` and `"TA_200m"` | Direct map via `CabauwAdapter` |
| **`:ameriflux`** | Contains `"TIMESTAMP_START"` or `"TA_H_V_R"` | Dynamic `stations.json` index lookup |
| **`:neon`** | Matches Regex `<prefix>_z\d+` (e.g., `temp_z1`) | Variable level extraction via `NEONAdapter` |
| **`:icos`** | Sparse 2-level headers (e.g., `T_10m`, `T_50m`) | MOST profile reconstruction via `ICOSAdapter` |
| **`:smear`** | Keys matching `HYY_` or `VAR_` prefixes | Direct mapping via `vars_lookup_compact.json` |

---

```julia


```

---

**Cross-Network Execution Example**

```julia
using DataFrames, Dates
using .UnifiedBLIngestion
using .UltraStability

# Example 1: NEON unstacked DataFrame
df_neon = DataFrame(
    datetime = [DateTime(2026, 8, 22, 12, 0)],
    temp_z1 = [22.1], temp_z2 = [21.5], temp_z3 = [20.8],
    ustar = [0.28], L = [-45.0]
)

# Ingest and extract Chebyshev coefficients automatically
neon_obs = ingest_boundary_layer_data(df_neon; prefix="temp", heights=[2.0, 10.0, 30.0], campaign="NEON-HARV")
neon_fingerprints = batch_spectral_fingerprints(neon_obs; n_coeffs=3)

# Example 2: ICOS 2-level sparse DataFrame (auto-upscaled via MOST)
df_icos = DataFrame(
    datetime = [DateTime(2026, 8, 22, 12, 0)],
    T_10m = [15.4], T_50m = [13.2],
    ustar = [0.18], L = [120.0]
)

icos_obs = ingest_boundary_layer_data(df_icos; campaign="ICOS-SE-S2")
icos_fingerprints = batch_spectral_fingerprints(icos_obs; n_coeffs=4)

```