# SBLToolkit.jl

A modular Julia pipeline for ingesting stable boundary layer (SBL) observational datasets, automated artifact cleansing, SVD modal signal isolation, and Geometric Surface Pattern Theory (GSPT) regime-transition inversions.

Designed to bridge field campaign NetCDF records with reduced-order atmospheric modeling in `GeoABL.jl` and WSINDy dynamical discovery pipelines.

---

## Key Capabilities

* **Catalog Ingestion**: Dynamic integration with `SpectralBL-Analytics/data/datasets.json` catalog specifications, auto-routing, missing-safe `coalesce` NetCDF parsing, and defensive `ispath` unmounted directory guards.
* **Artifact Cleansing & Defensive Ingestion**: Resolves $u_*$ flatline defaults ($0.3\text{ m/s}$), handles missing `NaN` entries without global stencil corruption via sub-grid derivative isolation, and applies Planar Fit coordinate rotations for snowpack/sloped terrain.
* **GSPT Transition Surface Mapping**: Evaluates the Phase 2 dynamic transition coordinate metric $R_{\text{coord}}(z, t)$ across space and time, identifying turbulence collapse boundaries ($Ri \approx Ri_{\text{cr}} \approx 0.25$).
* **Track A Primitive Regularization**: Applies Tikhonov-Morozov regularization to primitive fields ($u, v, \theta_v$) under a Modified Morozov Discrepancy Principle (MDP) to eliminate non-linear quotient singularities near Low-Level Jet (LLJ) noses.
* **Triple-Point Dispersion ($\delta_{\text{TP}}$) Quality Control**: Tracks non-dimensional spatial spread ($\delta_{\text{TP}}$) between diffusivity cutoffs ($z_K$), TKE level-sets ($z_e$), and TKE gradient extrema ($z_{e_z}$) to enforce automated quality-control gating ($\delta_{\text{TP}} < 0.10$).
* **Opt-In Diagnostic Synthesis**: Reconstructs sparse turbulence profiles via interior mast linear interpolation and standard boundary layer flux decay ($F(z) = F_0 \max\left(0, 1 - z/h_{\text{sbl}}\right)$) without contaminating baseline production runs.
* **Modal Decomposition**: Extracts SVD spatial profiles ($V(z)$) and temporal principal components ($PC_k(t)$) from velocity anomaly matrices $U(z, t) - \bar{U}(z)$.
* **Wave/Jet Separation**: Zero-phase Butterworth lowpass filtering isolates Low-Level Jets ($PC_{2,\text{LLJ}}$) from Internal Gravity Wave perturbations ($PC_{2,\text{IGW}}$) without boundary-edge artifacts.
* **Sub-Grid Core Tracking**: Recovers continuous jet-core height $z_{\text{LLJ}}(t)$ using 3-point parabolic peak interpolation.
* **Spectral Profile Encoding**: Multi-level profiles pass through logarithmic Chebyshev projection to generate compact spectral fingerprints ($c_0, c_1, c_2, c_3$).

---

## Ingested Campaign Benchmarks

| Campaign ID | Archetype | Specific Adjustments / Transforms | Diagnostic & Inversion Status |
| --- | --- | --- | --- |
| `floss_ii` | Marginal Stability | Planar Fit sonic rotation, $u_* = 0.3\text{ m/s}$ replacement | Guarded (`ispath` check cleanly skips if unmounted) |
| `bllast` | Convective-to-Stable | Level-specific temperature offsets, Flag-H/LE stationarity masking | Guarded (`ispath` check cleanly skips if unmounted) |
| `cases_99` | Brittle Transition | Standard double rotation, high-stability nocturnal profile extraction | **Full Inversion Success** ($N_z=6$, LLJ nose fold identified at $45\text{ m}$) |
| `gabls3` | Diurnal Process Benchmark | 24-hour diurnal cycle profile and jet tracking | **Diagnostic Synthesis Success** ($N_z=38$, 900 finite $R_{\text{coord}}$ points recovered) |
| `sheba` | Polar Snow/Ice SBL | Surface energy balance auditing, extreme stability Monin-Obukhov scaling | **Bypassed Guard** ($N_z=2 < 3$, finite-difference stencil under-determination) |

---

## Repository Architecture Map

| Subsystem / Path | Core Components | Primary Responsibilities |
| --- | --- | --- |
| **`src/`** | `GSPTTimeLoop.jl`, `GSPTPhase2.jl`, `gspt_triple_point.jl`, `gspt_nonuniform_gradient.jl`, `SBLToolkit.jl` | Houses GSPT Phase 2 inversion engine, Track A Tikhonov regularization, non-uniform derivative operators, and triple-point dispersion tracking. |
| **`src/ultra/adapters/`** | `ameriflux_adapter.jl`, `cabauw_adapter.jl`, `icos_adapter.jl`, `neon_adapter.jl`, `sheba_adapter.jl`, `bllast_adapter.jl` | Normalizes heterogeneous tower network files into standardized observation contracts. |
| **`src/ultra/`** | `core_types.jl`, `forcing.jl`, `stability.jl`, `spectral_engine.jl`, `spectral_reconstruction.jl`, `UnifiedBLIngestion.jl` | Houses core math engines for Chebyshev profile projection, SCM forcing generation, stability regimes, and auto-routing. |
| **`data/`** | `raw/` (`bllast`, `cases99`, `floss`, `gabls3`, `sheba`), `processed/`, `catalog.json` | Stores raw field campaign observations, transformed JLD2 datasets, and network site lookup registries. |
| **`scripts/`** | `plot_gspt_transition.jl`, `process_campaign.jl`, `convert_ncar_cases99.jl`, `convert_sheba.jl` | Operational CLI entrypoints for dataset transformations, dynamic heatmap generation, and stability audits. |
| **`reports/generated/`** | `gspt_phase2/` (`*_transition_heatmap.png`, `gspt_triple_point_dispersion.png`) | Generated diagnostic outputs, transition surface plots, and dispersion analysis artifacts. |
| **`notes/`** | `Stability Closure Comparison.md`, `audit_seb_nc.md`, Level 5 Scientific Validation Reports | Scientific audit logs, surface energy balance diagnostics, and regularization closure proofs. |

---

## Directory Structure

```text
SBLToolkit.jl/
├── Project.toml
├── Makefile
├── README.md
├── run_pipeline.sh
├── src/
│   ├── SBLToolkit.jl                 # Primary module interface
│   ├── Catalog.jl                    # Catalog JSON parser
│   ├── Transformations.jl            # Planar Fit & rotation math
│   ├── IngestBLLAST.jl               # Campaign-specific parsers
│   ├── GSPTPhase2.jl                 # GSPT coordinate metric inversion & Z0HR schemes
│   ├── GSPTTimeLoop.jl               # Defensive ingestion, auditing & spatial derivative loops
│   ├── gspt_triple_point.jl          # Non-dimensional triple-point dispersion pipeline
│   ├── gspt_nonuniform_gradient.jl  # Sub-grid finite difference operator builders
│   └── ultra/                        # Ultra-stability ingestion and analysis
│       ├── core_types.jl
│       ├── forcing.jl
│       ├── stability.jl
│       ├── spectral_engine.jl
│       ├── spectral_reconstruction.jl
│       ├── UnifiedBLIngestion.jl
│       └── adapters/                 # Network-specific adapters
├── scripts/                          # CLI transformation and analysis tools
│   └── plot_gspt_transition.jl       # Dynamic GSPT transition heatmap generator
├── reports/                          # Generated engineering reports and plots
│   └── generated/
│       └── gspt_phase2/
├── test/                             # Integration & diagnostic test suite
│   ├── runtests.jl
│   └── test_triple_point_qc.jl       # Automated δ_TP < 0.10 QC gating suite
└── data/
    ├── raw/                          # Input NetCDF/CSV files
    └── processed/                    # Output JLD2 matrices

```

---

## Pipeline Execution Workflow

1. **Data Standardization**: Conversion scripts in `scripts/convert_*.jl` read raw campaign files in `data/raw/` and feed them through `UnifiedBLIngestion.jl` using network adapters.
2. **GSPT Transition Heatmap Generation**: Run the Phase 2 inversion pipeline to generate spatial transition heatmaps ($R_{\text{coord}}(z, t)$):
```bash
# Run baseline production heatmaps
make gspt

# Run with full diagnostic logging and raw fallback plots
GSPT_DEBUG_AUDIT=1 GSPT_WRITE_RAW_DIAGNOSTIC=1 GSPT_SYNTHESIZE_MISSING_FLUXES=1 make gspt

```


3. **Triple-Point Dispersion Analysis**: Quantify physical turbulence collapse and verify spatial marker alignment ($\delta_{\text{TP}} \to 0$):
```bash
julia --project=. src/gspt_triple_point.jl

```


4. **Spectral & Modal Decomposition**: Multi-level profiles pass through `spectral_engine.jl` to compute logarithmic Chebyshev coefficients ($c_0, c_1, c_2, c_3$) and isolate velocity principal components.

---

## Verification & Testing

Execute the complete automated test suite (including sub-grid gradient operators and the $\delta_{\text{TP}} < 0.10$ quality-control gate):

```bash
julia --project=. test/runtests.jl

```