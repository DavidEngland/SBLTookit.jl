# SBLToolkit.jl

A modular Julia pipeline for ingesting stable boundary layer (SBL) observational datasets, automated artifact cleansing, and SVD modal signal isolation.

Designed to bridge field campaign NetCDF records with reduced-order atmospheric modeling in `GeoABL.jl` and WSINDy dynamical discovery pipelines.

---

## Key Capabilities

* **Catalog Ingestion**: Dynamic integration with `SpectralBL-Analytics/data/datasets.json` catalog specifications and auto-routing.
* **Artifact Cleansing**: Resolves $u_*$ flatline defaults (0.3 m/s), handles missing NaNs via row-wise linear interpolation, and applies Planar Fit coordinate rotations for snowpack/sloped terrain.
* **Modal Decomposition**: Extracts SVD spatial profiles ($V(z)$) and temporal principal components ($PC_k(t)$) from velocity anomaly matrices $U(z,t) - \bar{U}(z)$.
* **Wave/Jet Separation**: Zero-phase Butterworth lowpass filtering isolates Low-Level Jets ($PC_{2,\text{LLJ}}$) from Internal Gravity Wave perturbations ($PC_{2,\text{IGW}}$) without boundary-edge artifacts.
* **Sub-Grid Core Tracking**: Recovers continuous jet-core height $z_{\text{LLJ}}(t)$ using 3-point parabolic peak interpolation.
* **Spectral Profile Encoding**: Multi-level profiles pass through logarithmic Chebyshev projection to generate compact spectral fingerprints ($c_0, c_1, c_2, c_3$).

---

## Ingested Campaign Benchmarks

| Campaign ID | Archetype | Specific Adjustments / Transforms |
| --- | --- | --- |
| `floss_ii` | Marginal Stability | Planar Fit sonic rotation, $u_* = 0.3\text{ m/s}$ replacement |
| `bllast` | Convective-to-Stable | Level-specific temperature offsets, Flag-H/LE stationarity masking |
| `cases_99` | Brittle Transition | Standard double rotation, high-stability nocturnal profile extraction |
| `gabls3` | Diurnal Process Benchmark | 24-hour diurnal cycle profile and jet tracking |
| `sheba` | Polar Snow/Ice SBL | Surface energy balance auditing, extreme stability Monin-Obukhov scaling |

---

## Repository Architecture Map

| Subsystem / Path | Core Components | Primary Responsibilities |
| --- | --- | --- |
| **`src/ultra/adapters/`** | `ameriflux_adapter.jl`, `cabauw_adapter.jl`, `icos_adapter.jl`, `neon_adapter.jl`, `smear_adapter.jl`, `sheba_adapter.jl`, `bllast_adapter.jl` | Normalizes heterogeneous tower network files into standardized observation contracts. |
| **`src/ultra/`** | `core_types.jl`, `forcing.jl`, `stability.jl`, `spectral_engine.jl`, `spectral_reconstruction.jl`, `UnifiedBLIngestion.jl` | Houses core math engines for Chebyshev profile projection, SCM forcing generation, stability regimes, and auto-routing. |
| **`data/`** | `raw/` (`bllast`, `cases99`, `floss`, `gabls3`, `sheba`), `processed/`, `catalog.json` | Stores raw field campaign observations, transformed datasets, and network site lookup registries. |
| **`scripts/`** | `process_campaign.jl`, `convert_ncar_cases99.jl`, `convert_sheba.jl`, `RiCurvature.jl`, `etas.jl` | Operational CLI entrypoints for dataset transformations, stability audits, and bulk calculations. |
| **`notes/`** | `Stability Closure Comparison.md`, `audit_seb_nc.md`, `21Aug2026.md` | Scientific audit logs, surface energy balance diagnostics, and closure theory comparisons. |

---

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
│   ├── IngestBLLAST.jl       # Campaign-specific parsers
│   └── ultra/                # Ultra-stability ingestion and analysis
│       ├── core_types.jl
│       ├── forcing.jl
│       ├── stability.jl
│       ├── spectral_engine.jl
│       ├── spectral_reconstruction.jl
│       ├── UnifiedBLIngestion.jl
│       └── adapters/         # Network-specific adapters
├── scripts/                  # CLI transformation and analysis tools
├── test/                     # Integration test suite
├── notes/                    # Scientific audits and closure notes
└── data/
    ├── raw/                  # Input NetCDF/CSV files
    └── processed/            # Output JLD2 matrices

```

---

## Pipeline Execution Workflow

1. **Data Standardization**: Conversion scripts in `scripts/convert_*.jl` read raw campaign files in `data/raw/` and feed them through `UnifiedBLIngestion.jl` using network adapters.
2. **Spectral Analysis**: Multi-level profiles pass through `spectral_engine.jl` to compute logarithmic Chebyshev coefficients ($c_0, c_1, c_2, c_3$).
3. **SCM Driving**: `forcing.jl` generates Single Column Model boundary conditions (`theta_tendency`, `zeta_reference`, heat fluxes) categorized by regime using `stability.jl`.
4. **Modal Signal Isolation**: Velocity matrices undergo SVD decomposition to isolate principal components and track continuous jet-core dynamics.

---

## Verification & Testing

Execute the automated integration test suite:

```bash
julia --project=. test/runtests.jl

```