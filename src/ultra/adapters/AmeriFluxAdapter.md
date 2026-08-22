`AmeriFluxAdapter.jl` completes the observational ingestion pipeline by extending boundary layer standardization across the 60+ tower AmeriFlux network, handling site-specific sensor configurations via external JSON metadata.

**Multi-Site Ingestion Architecture Comparison**

| Adapter | Network | Sensor Heights | Metadata Source | Missing Value Sentinel |
| --- | --- | --- | --- | --- |
| **`CabauwAdapter`** | CESAR (Netherlands) | Static 7-node mast ($2\text{m} \to 200\text{m}$) | Hardcoded (`CabauwTower`) | `missing` / `NaN` |
| **`SmearAdapter`** | SMEAR (Finland) | Variable lookup keys ($4\text{m} \to 125\text{m}$) | `vars_lookup_compact.json` | `NaN` |
| **`AmeriFluxAdapter`** | AmeriFlux (North America) | Dynamic per-site ($H_1, H_2, \dots$) | `stations.json` | `-9999` |

---

**JSON Station Registry Schema (`data/ameriflux/stations.json`)**

To support dynamic height mapping for sites like `US-ARM` (Southern Great Plains) or `US-xCP`, the registry maps FP-standard indices ($1, 2, \dots, N$) to physical heights:

```json
{
  "sites": [
    {
      "id": "US-ARM",
      "name": "ARM Southern Great Plains site",
      "measurement_heights_m": [2.0, 10.0, 25.0, 60.0],
      "z0m": 0.05,
      "d_displacement": 0.0
    },
    {
      "id": "US-Ha1",
      "name": "Harvard Forest EMS Tower",
      "measurement_heights_m": [12.0, 18.0, 24.0, 29.0],
      "z0m": 1.10,
      "d_displacement": 15.2
    }
  ]
}

```

---

**Code Quality & Safety Refinement**

Two minor edge cases in `extract_ameriflux_observations` should be guarded against:

1. **Empty Bounds Indexing**: If a row contains `ZL` but all temperature levels are missing or filtered out, accessing `heights_m[1]` throws a `BoundsError`.
2. **Constructor Consistency**: Ensure `StandardizedBLObservation` constructor calls match the field order defined across `CabauwAdapter`.

```julia
# 1. Calculate L_obukhov ONLY after confirming valid heights exist
if n_valid >= min_levels
    L_obukhov = 9999.0
    if haskey(row, "ZL") && !ismissing(row["ZL"])
        zl = to_float(row["ZL"], NaN)
        if !isnan(zl) && abs(zl) > 1e-5
            # Safely use the first valid physical height
            L_calc = heights_m[1] / zl
            L_obukhov = clamp(L_calc, -9999.0, 9999.0)
        end
    end

    # 2. Instantiate standardized observation struct
    obs = StandardizedBLObservation(
        dt,
        "AMERIFLUX",
        heights_m,
        temp_values,
        ustar,
        L_obukhov,
        z0m,
        (n_valid >= 3),
        n_valid
    )
    push!(observations, obs)
end

```

Would you like to build a unified unified loader function (`load_campaign_observation`) that routes any site CSV/DataFrame (Cabauw, SMEAR, or AmeriFlux) into `SpectralEngine` automatically?