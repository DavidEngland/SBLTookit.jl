module SBLGating  
  
using LinearAlgebra  
  
export BifurcationGatingParams, apply_geometry_gating  
  
"""  
    BifurcationGatingParams{T}  
  
Parameters defining the fast-eigenvalue gating threshold and minimum residual mixing limits.  
"""  
struct BifurcationGatingParams{T<:AbstractFloat}  
    epsilon_lambda::T  # Eigenvalue threshold (e.g., 1e-4)  
    km_floor::T        # Floor eddy diffusivity for momentum (m^2/s)  
    kh_floor::T        # Floor eddy diffusivity for heat (m^2/s)  
end  
  
function BifurcationGatingParams(; epsilon_lambda=1e-4, km_floor=0.1, kh_floor=0.01)  
    return BifurcationGatingParams(epsilon_lambda, km_floor, kh_floor)  
end  
  
"""  
    apply_geometry_gating(lambda_f, K_m, K_h, params)  
  
Evaluates the fast eigenvalue condition λ_f ≥ -**ϵ**_λ. If satisfied, overrides standard  
Monin-Obukhov/diagnostic turbulence shutdowns to prevent false runaway decoupling (Pathology B).  
"""  
function apply_geometry_gating(  
    lambda_f::T,  
    K_m::T,  
    K_h::T,  
    params::BifurcationGatingParams{T}  
) where {T<:AbstractFloat}  
  
    # Flow remains dynamically active if fast eigenvalue is non-negative within tolerance  
    is_hyperbolic_active = lambda_f >= -params.epsilon_lambda  
  
    if is_hyperbolic_active  
        # Override diagnostic mixing shutdown; retain active mixing layer  
        K_m_gated = max(K_m, params.km_floor)  
        K_h_gated = max(K_h, params.kh_floor)  
        is_gated = true  
    else  
        # Permit physical turbulence decay/decoupling  
        K_m_gated = K_m  
        K_h_gated = K_h  
        is_gated = false  
    end  
  
    return K_m_gated, K_h_gated, is_gated  
end  
  
end # module SBLGating  
