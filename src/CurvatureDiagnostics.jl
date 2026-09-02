module CurvatureDiagnostics

using SmoothingSplines
using LinearAlgebra
using Statistics

export process_dns_profile, CurvatureOutput

struct CurvatureOutput
    z::Vector{Float64}
    Rig::Vector{Float64}
    chi::Vector{Float64}
    zeta::Vector{Float64}
    zeta_z::Vector{Float64}
    zeta_zz::Vector{Float64}
    C_M::Vector{Float64}
    classification::Vector{Symbol}
end

function spline_derivative(spline::SmoothingSpline{Float64}, x::Vector{Float64}, order::Int)
    order in (1, 2) || throw(ArgumentError("Derivative order must be 1 or 2."))

    derivative = Vector{Float64}(undef, length(x))
    knots = spline.Xdesign
    values = spline.g
    second_derivatives = spline.γ
    last_knot = length(knots)

    for (index, x_value) in pairs(x)
        left = clamp(searchsortedlast(knots, x_value), 1, last_knot - 1)
        right = left + 1
        interval = knots[right] - knots[left]
        t = clamp((x_value - knots[left]) / interval, 0.0, 1.0)
        left_second = left == 1 ? 0.0 : second_derivatives[left - 1]
        right_second = right == last_knot ? 0.0 : second_derivatives[right - 1]

        if order == 1
            derivative[index] = (values[right] - values[left]) / interval +
                                interval / 6 * ((1 - 3 * (1 - t)^2) * left_second +
                                                (3 * t^2 - 1) * right_second)
        else
            derivative[index] = (1 - t) * left_second + t * right_second
        end
    end

    return derivative
end

"""
    process_dns_profile(z, u_raw, v_raw, theta_raw; g=9.81, theta0=300.0, kappa=0.4,
                        tau_raw=nothing, heat_flux_raw=nothing)

Ingests raw DNS vertical profiles, fits GCV smoothing splines to primitive fields,
analytically evaluates spatial derivatives, and computes the C_M mapping fraction.
"""
function process_dns_profile(
    z::Vector{Float64},
    u_raw::Vector{Float64},
    v_raw::Vector{Float64},
    theta_raw::Vector{Float64};
    g::Float64 = 9.81,
    theta0::Float64 = 300.0,
    kappa::Float64 = 0.4,
    delta_fold::Float64 = 0.10,
    C_M_min::Float64 = 0.70,
    epsilon_C::Float64 = 10.0,
    smoothing_parameter::Float64 = 1.0,
    curvature_smoothing_parameter::Float64 = 1.0,
    tau_raw::Union{Nothing,AbstractVector{<:Real}} = nothing,
    heat_flux_raw::Union{Nothing,AbstractVector{<:Real}} = nothing,
    state_singularity::Union{Nothing,AbstractVector{<:Real}} = nothing,
    verified_saddle_node::Union{Nothing,AbstractVector{Bool}} = nothing,
    closure_residual::Union{Nothing,AbstractVector{<:Real}} = nothing,
    delta_obs::Float64 = Inf,
    state_singularity_tol::Float64 = 1e-8
)
    Nz = length(z)
    L_z = z[end] - z[1]

    all(length(values) == Nz for values in (u_raw, v_raw, theta_raw)) ||
        throw(ArgumentError("All profile vectors must have the same length as z."))
    (tau_raw === nothing) == (heat_flux_raw === nothing) ||
        throw(ArgumentError("tau_raw and heat_flux_raw must be supplied together."))
    delta_fold >= 0.0 || throw(ArgumentError("delta_fold must be non-negative."))
    0.0 <= C_M_min <= 1.0 || throw(ArgumentError("C_M_min must lie in [0, 1]."))
    delta_obs >= 0.0 || throw(ArgumentError("delta_obs must be non-negative."))
    state_singularity_tol >= 0.0 || throw(ArgumentError("state_singularity_tol must be non-negative."))

    for (name, values) in (("tau_raw", tau_raw),
                           ("heat_flux_raw", heat_flux_raw),
                           ("state_singularity", state_singularity),
                           ("verified_saddle_node", verified_saddle_node),
                           ("closure_residual", closure_residual))
        values === nothing || length(values) == Nz ||
            throw(ArgumentError("$name must have the same length as z."))
    end

    # Step 1: Primitive Field GCV Spline Smoothing (Enforcing Operator Non-Commutation)
    # Fit splines to raw primitive fields directly to avoid 1/U_z singularities
    spl_u = fit(SmoothingSpline, z, u_raw, smoothing_parameter)
    spl_v = fit(SmoothingSpline, z, v_raw, smoothing_parameter)
    spl_theta = fit(SmoothingSpline, z, theta_raw, smoothing_parameter)

    # Evaluate smoothed primitives and analytical derivatives
    u = predict(spl_u, z)
    u_z = spline_derivative(spl_u, z, 1)
    u_zz = spline_derivative(spl_u, z, 2)

    v = predict(spl_v, z)
    v_z = spline_derivative(spl_v, z, 1)
    v_zz = spline_derivative(spl_v, z, 2)

    theta = predict(spl_theta, z)
    theta_z = spline_derivative(spl_theta, z, 1)
    theta_zz = spline_derivative(spl_theta, z, 2)

    # Step 2: Construct Gradient Richardson Number Profile
    S2 = @. u_z^2 + v_z^2
    N2 = @. (g / theta0) * theta_z
    Rig = @. N2 / max(S2, 1e-12)

    # Step 3: Construct the inverse Obukhov profile before forming \zeta.
    # chi remains finite as L approaches +/-infinity; no L inversion is performed.
    chi = if tau_raw === nothing
        @. (g / theta0) * theta_z / max(S2, 1e-12)
    else
        tau = Float64.(tau_raw)
        heat_flux = Float64.(heat_flux_raw)
        all(isfinite, tau) && all(isfinite, heat_flux) ||
            throw(ArgumentError("tau_raw and heat_flux_raw must contain finite values."))
        all(tau .>= 0.0) || throw(ArgumentError("tau_raw must be non-negative."))
        @. -kappa * (g / theta0) * heat_flux / max(tau, eps(Float64))^(3 / 2)
    end
    zeta = @. z * chi

    # Differentiate the smoothed inverse profile, then apply the product rule.
    spl_chi = fit(SmoothingSpline, z, chi, curvature_smoothing_parameter)
    chi_z = spline_derivative(spl_chi, z, 1)
    chi_zz = spline_derivative(spl_chi, z, 2)
    zeta_z = @. chi + z * chi_z
    zeta_zz = @. 2.0 * chi_z + z * chi_zz

    # Step 4: Evaluate Closure Derivatives R'(\zeta) and R''(\zeta)
    # Using a smooth log-linear baseline closure R(\zeta) = \zeta * (1 + 5\zeta) / (1 + 15\zeta)^2
    R_prime = zeros(Nz)
    R_double_prime = zeros(Nz)
    for i in 1:Nz
        z_i = max(zeta[i], 0.0)
        # Analytical derivatives of standard smooth stability closure R(\zeta)
        R_prime[i] = (1.0 + 5.0 * z_i) / (1.0 + 15.0 * z_i)^2
        R_double_prime[i] = -25.0 * z_i / (1.0 + 15.0 * z_i)^3
    end

    # Step 5: Evaluate Local Adaptive Curvature Scale K_0 and Mapping Fraction C_M
    # Compute physical profile second derivative Rig_zz
    spl_Rig = fit(SmoothingSpline, z, Rig, curvature_smoothing_parameter)
    Rig_zz = spline_derivative(spl_Rig, z, 2)

    K0 = median(abs.(Rig_zz)) + 1e-6

    C_constitutive = @. abs(R_double_prime * (zeta_z^2))
    C_mapping = @. abs(R_prime * zeta_zz)
    C_M = @. C_mapping / (C_constitutive + C_mapping + epsilon_C * K0)

    # Step 6: Five-state classification matrix.
    # Missing state-space/closure evidence is treated as unverified, not as proof.
    state_singularity_values = state_singularity === nothing ? fill(Inf, Nz) : Float64.(state_singularity)
    verified_values = verified_saddle_node === nothing ? falses(Nz) : verified_saddle_node
    closure_values = closure_residual === nothing ? fill(0.0, Nz) : abs.(Float64.(closure_residual))

    classification = Vector{Symbol}(undef, Nz)
    for i in 1:Nz
        fold_scale = L_z * abs(zeta_zz[i])
        is_fold_prox = fold_scale > 0.0 && abs(zeta_z[i]) <= delta_fold * fold_scale
        is_mapping_dom = C_M[i] >= C_M_min
        is_singular = state_singularity_values[i] <= state_singularity_tol
        is_verified = verified_values[i]
        is_closure_valid = closure_values[i] <= delta_obs

        if !is_closure_valid
            classification[i] = :ClosureBreakdown
        elseif is_verified && is_fold_prox
            classification[i] = :HybridFold
        elseif is_verified
            classification[i] = :VerifiedSaddleNode
        elseif is_fold_prox && is_mapping_dom && !is_singular
            classification[i] = :PureCoordinateFold
        else
            classification[i] = :DynamicSingularityCandidate
        end
    end

    return CurvatureOutput(z, Rig, chi, zeta, zeta_z, zeta_zz, C_M, classification)
end

end # module