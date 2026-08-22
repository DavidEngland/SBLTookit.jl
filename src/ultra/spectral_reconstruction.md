The inverse spectral transformation maps Chebyshev coefficients $c_k$ back onto a dense 3D SCM vertical grid $z_{\text{grid}}$ by evaluating the synthesis matrix $V_{\text{target}}$ at the normalized layer target heights $x_{\text{grid}} \in [-1, 1]$.


**End-to-End Tower-to-SCM Interpolation Pipeline**

```julia
using .CoreTypes
using .SpectralEngine
using .SpectralReconstruction

# 1. Sparse multi-height tower data from SMEAR (e.g. Hyytiälä 7-level temperature)
tower_z = [4.2, 8.4, 16.8, 33.6, 50.4, 67.2, 125.0]
tower_T = [15.2, 14.8, 14.3, 13.9, 13.5, 13.2, 12.1]
d_0 = 12.0 # Pine forest displacement height

meta = ProfileMetadata(d_0)
prof = MeteorologicalProfile(tower_z, tower_T, meta)

# 2. Compress profile into 4 spectral modes (forward projection)
c_modes = chebyshev_fingerprint(prof; n_coeffs=4, height_mapping=:log)

# 3. Reconstruct onto a dense 50-level SCM model grid across [13.0m, 125.0m]
scm_grid_heights = collect(range(13.0, 125.0, length=50))
tower_bounds = (extrema(tower_z)...,)

scm_temperature_profile = reconstruct_from_chebyshev(
    c_modes,
    scm_grid_heights,
    tower_bounds;
    height_mapping=:log,
    canopy_displacement=d_0
)

```

**Key Advantages**

* **Smooth Derivatives:** Prevents standard linear or spline interpolation oscillations near boundary layer transitions.
* **Mass Conservation:** Preserves integral layer averages ($c_0$) while interpolating between non-uniform tower sensors.