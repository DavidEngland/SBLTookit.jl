#!/usr/bin/env julia
#SBLTookit.jl scripts/heatmaps.jl
using Plots
using ColorSchemes

# -------------------------------------------------------------------
# 1. Non-Uniform 1D Finite Difference Derivative Stencil
# -------------------------------------------------------------------
"""
    non_uniform_gradient(f::Vector{Float64}, z::Vector{Float64})

Computes 2nd-order accurate 1D spatial derivative df/dz on a non-uniformly
spaced vertical grid `z`.
"""
function non_uniform_gradient(f::Vector{Float64}, z::Vector{Float64})
    n = length(z)
    df = zeros(Float64, n)

    # 1st-order forward/backward differences at grid endpoints
    df[1] = (f[2] - f[1]) / (z[2] - z[1])
    df[end] = (f[end] - f[end-1]) / (z[end] - z[end-1])

    # 2nd-order central difference scheme for non-uniform interior nodes
    for i in 2:(n - 1)
        h1 = z[i] - z[i-1]
        h2 = z[i+1] - z[i]
        df[i] = (f[i+1] * h1^2 - f[i-1] * h2^2 + f[i] * (h2^2 - h1^2)) / (h1 * h2 * (h1 + h2))
    end
    return df
end

# -------------------------------------------------------------------
# 2. Synthetic Multi-Level Tower Dataset (24-Hour Diurnal Cycle)
# -------------------------------------------------------------------
const g = 9.81         # Gravitational acceleration [m/s^2]
const k_vk = 0.40      # von Kármán constant
const θ0 = 300.0       # Reference potential temperature [K]

# Tower measurement heights [m] (clustered near surface, wider aloft)
z_tower = [1.5, 3.0, 6.0, 10.0, 20.0, 30.0, 50.0, 80.0, 120.0]
nz = length(z_tower)

# Time axis (24 hours at 10-minute resolution)
hours = range(0.0, 24.0, length=145)
nt = length(hours)

# Preallocate flux fields
u_star = zeros(nz, nt)
wθ_flux = zeros(nz, nt)

for (j, t) in enumerate(hours)
    is_daytime = (t >= 6.0 && t <= 18.0)

    for (i, z) in enumerate(z_tower)
        if is_daytime
            # Daytime convective plume layer with weak flux divergence
            solar_phase = sin(π * (t - 6.0) / 12.0)
            wθ_flux[i, j] = 0.25 * solar_phase * exp(-z / 200.0)
            u_star[i, j] = 0.35 + 0.10 * solar_phase
        else
            # Stable nocturnal boundary layer with pronounced vertical divergence
            wθ_flux[i, j] = -0.04 * exp(-z / 25.0)  # Decaying downward heat flux
            u_star[i, j] = 0.18 * exp(-z / 80.0)    # Friction velocity reduction aloft
        end
    end
end

# -------------------------------------------------------------------
# 3. Compute Local L(z,t), ζ(z,t), and Non-Uniform Gradient ∂ζ/∂z
# -------------------------------------------------------------------
L_local = zeros(nz, nt)
zeta = zeros(nz, nt)
zeta_z = zeros(nz, nt)

for j in 1:nt
    for i in 1:nz
        flux = wθ_flux[i, j]
        if abs(flux) < 1e-6
            L_local[i, j] = sign(flux) * 1e5
        else
            L_local[i, j] = - (u_star[i, j]^3 * θ0) / (k_vk * g * flux)
        end
        zeta[i, j] = z_tower[i] / L_local[i, j]
    end

    # Compute spatial gradient ∂ζ/∂z using non-uniform finite differences
    zeta_z[:, j] = non_uniform_gradient(zeta[:, j], z_tower)
end

# Symmetric log-transform to compress unstable and stable scales
zeta_transformed = @. sign(zeta) * log10(1.0 + abs(zeta))

# -------------------------------------------------------------------
# 4. Publication-Quality Hövemöller Plot Generation
# -------------------------------------------------------------------
plt = heatmap(
    hours, z_tower, zeta_transformed,
    clims=(-1.5, 1.5),
    color=:balance,
    xlabel="Time of Day [UTC / Hours]",
    ylabel="Measurement Height z [m]",
    title="Hövemöller Diagram: Local Stability Parameter sgn(ζ) log₁₀(1 + |ζ|)",
    titlefontsize=11,
    guidefontsize=10,
    tickfontsize=9,
    colorbar_title="  sgn(ζ) log₁₀(1 + |ζ|)",
    xlims=(0, 24),
    xticks=0:3:24,
    yticks=[1.5, 10, 20, 30, 50, 80, 120],
    framestyle=:box,
    size=(900, 480),
    dpi=300
)

# Contour overlay for key local stability boundaries
contour!(
    hours, z_tower, zeta,
    levels=[0.0, 0.1, 0.5, 1.0, 2.0],
    color=:black,
    linewidth=0.8,
    linestyle=:dash,
    colorbar=false
)

savefig("hoevemoeller_stability_strata.png")
display(plt)