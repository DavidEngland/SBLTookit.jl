using LinearAlgebra  
using SBLGating  
using GABLS3Adapters  
  
function run_gabls3_step!(  
    E::Vector{T},  
    S::Vector{T},  
    N2::Vector{T},  
    zeta_z::Vector{T},  
    zeta_zz::Vector{T},  
    K_m::Vector{T},  
    K_h::Vector{T},  
    gating_states::Vector{GatingState{T}},  
    params::BifurcationGatingParams{T},  
    config::GSPTModelConfig{T};  
    enable_gating::Bool = true  
) where {T<:AbstractFloat}  
  
    nz = length(E)  
    gated_count = 0  
  
    for z in 1:nz  
        # 1. Linearize local fast-slow system  
        J_fast = map_gabls3_to_jacobian(z, E, S, N2, config)  
  
        # 2. Extract fast mode using eigenvector continuation  
        lambda_f = extract_fast_eigenvalue(J_fast, gating_states[z])  
  
        # 3. Calculate Fold Ratio audit metric  
        fold_ratio = compute_fold_ratio(zeta_z[z], zeta_zz[z], T(0.01))  
  
        if enable_gating  
            # 4. Evaluate hysteresis & geometric gating condition  
            Km_gated, Kh_gated, is_gated = update_gating_state!(  
                gating_states[z],  
                lambda_f,  
                zeta_z[z],  
                K_m[z],  
                K_h[z],  
                params  
            )  
  
            # Secondary Audit: Confirm mapping dominance when gated  
            if is_gated  
                K_m[z] = Km_gated  
                K_h[z] = Kh_gated  
                gated_count += 1  
            end  
        end  
    end  
  
    return gated_count  
end  
