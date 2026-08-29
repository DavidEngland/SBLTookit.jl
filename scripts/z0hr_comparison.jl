# ==============================================================================
# Technical Simulation: Zero-Offset Hyperbolic Regularization (Z0HR)
# Comparison of Smooth NWP Closures vs. Piecewise Classical Closures
# ==============================================================================
#
# This script computes, visualizes, and compares the regularized Z0HR momentum 
# (S_m) and heat (S_h) stability functions against traditional, non-differentiable 
# piecewise closures as a function of the Richardson number (Ri).
#
# It generates a publication-quality 2x2 multi-panel figure displaying:
#   1. Stability Functions Sm(Ri) and Sh(Ri) for both schemes.
#   2. First Derivatives (dSm/dRi, dSh/dRi) highlighting derivative continuity.
#   3. Residuals (Smooth - Classical) showing the O(ε) regularization behavior.
#   4. Log-scale Asymptotic Decay under strong stability (Ri > Ri_c).
#
# Dependencies: Plots.jl (standard Julia plotting library)
# To run this script locally:
#   1. open julia
#   2. import Pkg; Pkg.add("Plots")
#   3. include("z0hr_comparison.jl")
# ==============================================================================

using Plots
using LinearAlgebra

# ------------------------------------------------------------------------------
# 1. PARAMETERS & CONFIGURATION
# ------------------------------------------------------------------------------
const Ri_c  = 0.20   # Critical Richardson number (Businger-Dyer style)
const α_θ   = 1.00   # Inverse of turbulent Prandtl number prefactor (Pr_0^-1)
const Bu_m  = 4.70   # Momentum unstable stability constant (field-calibrated)
const Bu_h  = 4.70   # Heat unstable stability constant (field-calibrated)
const ε     = 1e-3   # Hyperbolic smoothing parameter (NWP default)

# ------------------------------------------------------------------------------
# 2. ALGEBRAIC PRIMITIVES & COORDINATE TRANSFORMATIONS
# ------------------------------------------------------------------------------

# Smooth max/min algebraic primitives using hyperbolic regularization
smooth_max(x, eps_param) = 0.5 * (x + sqrt(x^2 + eps_param^2))
smooth_min(x, eps_param) = 0.5 * (x - sqrt(x^2 + eps_param^2))

# Classical exact max/min primitives
class_max(x) = max(0.0, x)
class_min(x) = min(0.0, x)

# Regularized Richardson coordinates
get_Ri_plus(Ri, eps_param)  = 0.5 * (Ri + sqrt(Ri^2 + eps_param^2))
get_Ri_minus(Ri, eps_param) = 0.5 * (Ri - sqrt(Ri^2 + eps_param^2))

# ------------------------------------------------------------------------------
# 3. CORE CLOSURE INSTANTIATIONS
# ------------------------------------------------------------------------------

# --- A. Traditional Piecewise Classical Closures ---
function S_momentum_classical(Ri)
    Ri_plus  = class_max(Ri)
    Ri_minus = class_min(Ri)
    
    # Linear stable factor with sharp zero-cutoff
    f_stable = class_max(1.0 - Ri_plus / Ri_c)
    
    # Quadratic stable clamp & square-root unstable scaling
    S_m = (f_stable^2) * sqrt(1.0 - Bu_m * Ri_minus)
    return S_m
end

function S_heat_classical(Ri)
    Ri_plus  = class_max(Ri)
    Ri_minus = class_min(Ri)
    
    # Linear stable factor with sharp zero-cutoff
    f_stable = class_max(1.0 - Ri_plus / Ri_c)
    
    # Quadratic stable clamp & 0.75-power unstable scaling
    S_h = (1.0 / α_θ) * (f_stable^2) * ((1.0 - Bu_h * Ri_minus)^0.75)
    return S_h
end

# --- B. Zero-Offset Hyperbolic Regularization (Z0HR) ---
function S_momentum_smooth(Ri, eps_param)
    Ri_plus  = get_Ri_plus(Ri, eps_param)
    Ri_minus = get_Ri_minus(Ri, eps_param)
    
    # C^∞ smooth stable factor using hyperbolic smooth_max
    raw_factor = 1.0 - Ri_plus / Ri_c
    f_stable   = smooth_max(raw_factor, eps_param)
    
    # Smooth momentum function
    S_m = (f_stable^2) * sqrt(1.0 - Bu_m * Ri_minus)
    return S_m
end

function S_heat_smooth(Ri, eps_param)
    Ri_plus  = get_Ri_plus(Ri, eps_param)
    Ri_minus = get_Ri_minus(Ri, eps_param)
    
    # C^∞ smooth stable factor using hyperbolic smooth_max
    raw_factor = 1.0 - Ri_plus / Ri_c
    f_stable   = smooth_max(raw_factor, eps_param)
    
    # Smooth heat function
    S_h = (1.0 / α_θ) * (f_stable^2) * ((1.0 - Bu_h * Ri_minus)^0.75)
    return S_h
end

# ------------------------------------------------------------------------------
# 4. NUMERICAL EVALUATION OVER DETAILED Ri DOMAIN
# ------------------------------------------------------------------------------
# Grid of Richardson numbers spanning convective (unstable) to strongly stable
const Ri_grid = range(-0.5, stop=0.4, length=1000)
const h_diff  = 1e-6 # Step size for central finite differences

# Pre-allocate arrays
S_m_class = zeros(length(Ri_grid))
S_h_class = zeros(length(Ri_grid))
S_m_smooth = zeros(length(Ri_grid))
S_h_smooth = zeros(length(Ri_grid))

# Compute functions
for i in eachindex(Ri_grid)
    Ri = Ri_grid[i]
    S_m_class[i] = S_momentum_classical(Ri)
    S_h_class[i] = S_heat_classical(Ri)
    S_m_smooth[i] = S_momentum_smooth(Ri, ε)
    S_h_smooth[i] = S_heat_smooth(Ri, ε)
end

# Compute Residuals (Smooth Regularized - Piecewise Classical)
res_m = S_m_smooth .- S_m_class
res_h = S_h_smooth .- S_h_class

# Compute numerical derivatives (dSm/dRi and dSh/dRi) via central differences
dSm_dRi_class  = zeros(length(Ri_grid))
dSh_dRi_class  = zeros(length(Ri_grid))
dSm_dRi_smooth = zeros(length(Ri_grid))
dSh_dRi_smooth = zeros(length(Ri_grid))

for i in eachindex(Ri_grid)
    Ri = Ri_grid[i]
    
    # Central difference for momentum
    dSm_dRi_class[i]  = (S_momentum_classical(Ri + h_diff) - S_momentum_classical(Ri - h_diff)) / (2.0 * h_diff)
    dSm_dRi_smooth[i] = (S_momentum_smooth(Ri + h_diff, ε) - S_momentum_smooth(Ri - h_diff, ε)) / (2.0 * h_diff)
    
    # Central difference for heat
    dSh_dRi_class[i]  = (S_heat_classical(Ri + h_diff) - S_heat_classical(Ri - h_diff)) / (2.0 * h_diff)
    dSh_dRi_smooth[i] = (S_heat_smooth(Ri + h_diff, ε) - S_heat_smooth(Ri - h_diff, ε)) / (2.0 * h_diff)
end

# ------------------------------------------------------------------------------
# 5. MULTI-PANEL VISUALIZATION ARCHITECTURE (Plots.jl)
# ------------------------------------------------------------------------------
println("Generating publication-quality visualization panels...")
gr() # Use GR backend for fast, high-quality rendering

# Theme configuration
plot_font = "Helvetica"
default(fontfamily=plot_font, legendfontsize=9, titlefontsize=11, tickfontsize=9)

# Panel A: Stability Functions Sm and Sh
p1 = plot(title="A. Stability Functions (Sm, Sh)", xlabel="Richardson Number (Ri)", ylabel="Stability Function Value")
plot!(p1, Ri_grid, S_m_class,  label="Sm Piecewise Classical", color=:blue,  linestyle=:dash, width=1.5)
plot!(p1, Ri_grid, S_m_smooth, label="Sm Smooth Z0HR",        color=:blue,  linestyle=:solid, width=2.0)
plot!(p1, Ri_grid, S_h_class,  label="Sh Piecewise Classical", color=:orange,linestyle=:dash, width=1.5)
plot!(p1, Ri_grid, S_h_smooth, label="Sh Smooth Z0HR",        color=:orange,linestyle=:solid, width=2.0)
vline!(p1, [0.0], label="Neutral (Ri=0)", color=:gray, linestyle=:dot, width=1.2)
vline!(p1, [Ri_c], label="Critical (Ri_c=0.2)", color=:red, linestyle=:dot, width=1.2)

# Panel B: First Derivatives (dSm/dRi, dSh/dRi) showing derivative smoothing
p2 = plot(title="B. First Derivatives (dS/dRi)", xlabel="Richardson Number (Ri)", ylabel="First Derivative")
plot!(p2, Ri_grid, dSm_dRi_class,  label="dSm/dRi Classical", color=:blue,  linestyle=:dash, width=1.5)
plot!(p2, Ri_grid, dSm_dRi_smooth, label="dSm/dRi Z0HR",      color=:blue,  linestyle=:solid, width=2.0)
plot!(p2, Ri_grid, dSh_dRi_class,  label="dSh/dRi Classical", color=:orange, linestyle=:dash, width=1.5)
plot!(p2, Ri_grid, dSh_dRi_smooth, label="dSh/dRi Z0HR",      color=:orange, linestyle=:solid, width=2.0)
vline!(p2, [0.0], label=false, color=:gray, linestyle=:dot, width=1.2)
vline!(p2, [Ri_c], label=false, color=:red, linestyle=:dot, width=1.2)

# Panel C: Residuals (Smooth Regularized - Piecewise Classical)
p3 = plot(title="C. Residual Analysis (Smooth - Classical)", xlabel="Richardson Number (Ri)", ylabel="Residual Value")
plot!(p3, Ri_grid, res_m, label="Sm Residual", color=:blue, linestyle=:solid, width=2.0)
plot!(p3, Ri_grid, res_h, label="Sh Residual", color=:orange, linestyle=:solid, width=2.0)
vline!(p3, [0.0], label=false, color=:gray, linestyle=:dot, width=1.2)
vline!(p3, [Ri_c], label=false, color=:red, linestyle=:dot, width=1.2)

# Panel D: Strong Stability Asymptotics (Ri > Ri_c) on Log Scale
stable_grid = range(0.15, stop=0.5, length=1000)
S_m_class_stable = [S_momentum_classical(r) for r in stable_grid]
S_m_smooth_stable = [S_momentum_smooth(r, ε) for r in stable_grid]
S_h_class_stable = [S_heat_classical(r) for r in stable_grid]
S_h_smooth_stable = [S_heat_smooth(r, ε) for r in stable_grid]

# Use 1e-12 as minimum floor to avoid log10(0) issues for classical
clamp_zero(x) = max(1e-12, x)

p4 = plot(title="D. Stable Asymptotics (Ri > Ri_c)", xlabel="Richardson Number (Ri)", ylabel="log10(Stability Function)")
plot!(p4, stable_grid, log10.(clamp_zero.(S_m_class_stable)),  label="Sm Classical (cutoff)", color=:blue,  linestyle=:dash, width=1.5)
plot!(p4, stable_grid, log10.(clamp_zero.(S_m_smooth_stable)), label="Sm Z0HR (asymptotic)",  color=:blue,  linestyle=:solid, width=2.0)
plot!(p4, stable_grid, log10.(clamp_zero.(S_h_class_stable)),  label="Sh Classical (cutoff)", color=:orange, linestyle=:dash, width=1.5)
plot!(p4, stable_grid, log10.(clamp_zero.(S_h_smooth_stable)), label="Sh Z0HR (asymptotic)",  color=:orange, linestyle=:solid, width=2.0)
vline!(p4, [Ri_c], label=false, color=:red, linestyle=:dot, width=1.2)

# Assemble multi-panel layout
total_plot = plot(p1, p2, p3, p4, layout=(2,2), size=(1000, 800))

# ------------------------------------------------------------------------------
# 6. OUTPUT & SUMMARY STATS
# ------------------------------------------------------------------------------
# Save the plot
output_png = "z0hr_stability_comparison.png"
png(total_plot, output_png)
println("Successfully saved comparison plot to: ", output_png)

# Compute exact numerical statistics to display
idx_neutral = findmin(abs.(Ri_grid))[2]
idx_critical = findmin(abs.(Ri_grid .- Ri_c))[2]

println("\n========================================================")
println("           Z0HR REGULARIZATION COMPARISON STATS")
println("========================================================")
println("Smoothing Parameter (ε): ", ε)
println("Critical Richardson Number (Ri_c): ", Ri_c)
println("--------------------------------------------------------")
println("Neutral Point (Ri = 0.0):")
println("  - Classical Sm:  ", S_m_class[idx_neutral])
println("  - Regularized Sm:", S_m_smooth[idx_neutral], " (Offset: ", res_m[idx_neutral], ")")
println("  - Classical Sh:  ", S_h_class[idx_neutral])
println("  - Regularized Sh:", S_h_smooth[idx_neutral], " (Offset: ", res_h[idx_neutral], ")")
println("--------------------------------------------------------")
println("Physical Cutoff Point (Ri = Ri_c):")
println("  - Classical Sm:  ", S_m_class[idx_critical])
println("  - Regularized Sm:", S_m_smooth[idx_critical], " (C^1 residual: ", res_m[idx_critical], ")")
println("  - Classical Sh:  ", S_h_class[idx_critical])
println("  - Regularized Sh:", S_h_smooth[idx_critical], " (C^1 residual: ", res_h[idx_critical], ")")
println("========================================================\n")
