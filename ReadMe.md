# SBLToolkit.jl

A modular Julia pipeline for ingesting stable boundary layer (SBL) observational datasets, automated artifact cleansing, and SVD modal signal isolation.

Designed to bridge field campaign NetCDF records with reduced-order atmospheric modeling in `GeoABL.jl` and WSINDy dynamical discovery pipelines.

## Key Capabilities

* **Catalog Ingestion**: Dynamic integration with `SpectralBL-Analytics/data/datasets.json` catalog specifications.
* **Artifact Cleansing**: Resolves $u_*$ flatline defaults (0.3 m/s), handles missing NaNs via row-wise linear interpolation, and applies Planar Fit coordinate rotations for snowpack/sloped terrain.
* **Modal Decomposition**: Extracts SVD spatial profiles ($V(z)$) and temporal principal components ($PC_k(t)$) from velocity anomaly matrices $U(z,t) - \bar{U}(z)$.
* **Wave/Jet Separation**: Zero-phase Butterworth lowpass filtering isolates Low-Level Jets ($PC_{2,\text{LLJ}}$) from Internal Gravity Wave perturbations ($PC_{2,\text{IGW}}$) without boundary-edge artifacts.
* **Sub-Grid Core Tracking**: Recovers continuous jet-core height $z_{\text{LLJ}}(t)$ using 3-point parabolic peak interpolation.

## Ingested Campaign Benchmarks

| Campaign ID | Archetype | Specific Adjustments / Transforms |
| :--- | :--- | :--- |
| `floss_ii` | Marginal Stability | Planar Fit sonic rotation, $u_* = 0.3\text{ m/s}$ replacement |
| `bllast` | Convective-to-Stable | Level-specific temperature offsets, Flag-H/LE stationarity masking |
| `cases_99` | Brittle Transition | Standard double rotation, high-stability nocturnal profile extraction |
| `gabls3` | Diurnal Process Benchmark | 24-hour diurnal cycle profile and jet tracking |

## Directory Structure

```text
SBLToolkit.jl/
├── Project.toml
├── Makefile
├── README.md
├── run_pipeline.sh
├── src/
│   ├── SBLToolkit.jl         # Primary module interface
│   ├── Catalog.jl            # Catalog JSON parser
│   ├── Transformations.jl    # Planar Fit & rotation math
│   └── IngestBLLAST.jl       # Campaign-specific parsers
├── scripts/
│   └── process_campaign.jl   # CLI worker for single dataset
└── data/
    ├── raw/                  # Input NetCDF files
    └── processed/            # Output JLD2 matrices