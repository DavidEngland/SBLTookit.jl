`SpectralEngine` projects non-uniformly sampled vertical tower profiles into a compact vector of orthogonal Chebyshev coefficients ($c_k$) using a pseudo-inverse Vandermonde matrix ($V^+$). Logarithmic height mapping ($x = 2 \frac{\ln z' - \ln z_{\min}}{\ln z_{\max} - \ln z_{\min}} - 1$) concentrates basis-function resolution near the surface where logarithmic gradients dominate.

**Chebyshev Coefficient Physical Interpretation**

| Coefficient | Polynomial $T_k(x)$ | Atmospheric Profile Feature |
| --- | --- | --- |
| $c_0$ | $1$ | Mean domain value across the tower layer |
| $c_1$ | $x$ | Primary vertical gradient (lapse rate / wind shear) |
| $c_2$ | $2x^2 - 1$ | Profile curvature (canopy top inversion / nocturnal decoupling) |
| $c_3$ | $4x^3 - 3x$ | Higher-order inflection (e.g., internal boundary layer transitions) |

---

**SMEAR Tower Profile Processing Example**

This snippet extracts temperature observations across 7 height levels at Hyytiälä (`station_id: 2`), constructs a `MeteorologicalProfile`, and generates a 4-mode Chebyshev fingerprint:

```julia
using .CoreTypes
using .SpectralEngine

# Mock definitions matching CoreTypes dependencies
struct ProfileMetadata
    canopy_displacement::Float64
end

struct MeteorologicalProfile
    heights::Vector{Float64}
    values::Vector{Float64}
    metadata::ProfileMetadata
end

"""
    extract_hyytiala_temperature_fingerprint(smear_frame::Dict{String, Float64})

Extracts Hyytiälä multi-height temperature variables, maps to displacement-adjusted
heights, and returns the log-spaced Chebyshev spectral coefficients.
"""
function extract_hyytiala_temperature_fingerprint(smear_frame::Dict{String, Float64})
    # Hyytiälä tower measurement heights (meters) and variable lookup keys
    tower_heights = [4.2, 8.4, 16.8, 33.6, 50.4, 67.2, 125.0]
    keys = ["HYY_META.T42", "HYY_META.T84", "HYY_META.T168",
            "HYY_META.T336", "HYY_META.T504", "HYY_META.T672", "HYY_META.T168"] # fallback

    # Reconstruct vertical profile array
    t_profile = [get(smear_frame, key, NaN) for key in keys]

    # Assume canopy zero-plane displacement d_0 = 12.0m for Hyytiälä pine forest
    meta = ProfileMetadata(12.0)

    # Filter valid heights strictly above displacement height d_0
    valid_idx = findall(h -> h > meta.canopy_displacement && !isnan(t_profile[h]), tower_heights)

    prof = MeteorologicalProfile(
        tower_heights[valid_idx],
        t_profile[valid_idx],
        meta
    )

    # Compute 4-coefficient fingerprint using log-height mapping
    return chebyshev_fingerprint(prof; n_coeffs=4, height_mapping=:log)
end

```

Would you like to build an inverse spectral reconstruction function to recreate 3D SCM vertical grid profiles from these Chebyshev coefficients?