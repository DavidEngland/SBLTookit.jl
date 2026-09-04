module GABLS3Adapters  
  
using LinearAlgebra  
using SBLGating  
  
export GSPTModelConfig, map_gabls3_to_jacobian, compute_fold_ratio  
  
"""  
    GSPTModelConfig{T}  
  
Configuration container for GSPT closure parameters and scale separation ratio.  
"""  
struct GSPTModelConfig{T<:AbstractFloat}  
    epsilon::T     # Timescale separation parameter (τ_fast / τ_slow ≈ 1e-3)  
    epsilon_c::T   # TKE regularization scale δ_reg (m²/s²)  
    l0::T          # Master turbulent length scale (m)  
end  
  
function GSPTModelConfig(; epsilon=1e-3, epsilon_c=1e-3, l0=15.0)  
    return GSPTModelConfig(epsilon, epsilon_c, l0)  
end  
  
"""  
    map_gabls3_to_jacobian(z_idx, E, S, N2, config)  
  
Constructs the 2x2 fast-slow Jacobian matrix at height index `z_idx`:  
    J = [ (1/ε) ∂F/∂e   (1/ε) ∂F/∂S ;  
               ∂g/∂e          ∂g/∂S ]  
"""  
function map_gabls3_to_jacobian(  
    z_idx::Int,  
    E::Vector{T},  
    S::Vector{T},  
    N2::Vector{T},  
    config::GSPTModelConfig{T}  
) where {T<:AbstractFloat}  
  
    # 1. State extraction with physical positivity floor  
    e_local = max(E[z_idx], T(1e-6))  
    s_local = max(S[z_idx], T(1e-6))  
    n2_local = N2[z_idx]  
  
    # 2. Hill-type buoyancy flux derivative: B(e) = B_max * e² / (e² + δ_reg²)  
    δ_reg = config.epsilon_c  
    B_max = config.l0 * n2_local  
    dB_de = B_max * (2 * e_local * δ_reg^2) / (e_local^2 + δ_reg^2)^2  
  
    # 3. Kolmogorov dissipation derivative: ε_diss = e^(3/2) / l0  
    dDiss_de = (1.5 * sqrt(e_local)) / config.l0  
  
    # 4. Fast vector field partials (∂F/∂e, ∂F/∂S)  
    Fe = config.l0 * s_local^2 - dB_de - dDiss_de  
    Fs = 2 * config.l0 * e_local * s_local  
  
    # 5. Slow vector field partials (∂g/∂e, ∂g/∂S)  
    gamma_s = T(0.1)  
    r_s = T(0.01)  
    ge = -gamma_s * s_local  
    gs = -gamma_s * e_local - r_s  
  
    # 6. Assemble scaled 2x2 Jacobian  
    inv_eps = one(T) / config.epsilon  
    return [ Fe * inv_eps   Fs * inv_eps ;  
             ge             gs           ]  
end  
  
"""  
    compute_fold_ratio(zeta_z, zeta_zz, C_constitutive)  
  
Calculates the Fold Ratio C_mapping = |curvature_mapping| / (|curvature_mapping| + |curvature_constitutive|).  
A value > 0.99 confirms that local profile curvature is overwhelmingly driven by coordinate mapping.  
"""  
function compute_fold_ratio(zeta_z::T, zeta_zz::T, C_constitutive::T) where {T<:AbstractFloat}  
    # Mapping curvature proxy scaled by gradient compression  
    C_mapping = abs(zeta_zz) / (abs(zeta_z) + T(1e-5))  
    total_curvature = C_mapping + abs(C_constitutive) + T(1e-10)  
    return C_mapping / total_curvature  
end  
  
end # module GABLS3Adapters  
