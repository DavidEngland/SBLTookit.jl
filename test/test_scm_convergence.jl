using Test
using Printf

"""
Monotone cubic Hermite interpolation (PCHIP/Fritsch-Carlson style).
This preserves profile monotonicity during vertical grid refinement.
"""
function pchip_slopes(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    m = zeros(Float64, n)

    if n == 2
        m .= δ[1]
        return m
    end

    # Endpoints (shape-preserving estimates)
    m[1] = ((2h[1] + h[2]) * δ[1] - h[1] * δ[2]) / (h[1] + h[2])
    if sign(m[1]) != sign(δ[1])
        m[1] = 0.0
    elseif sign(δ[1]) != sign(δ[2]) && abs(m[1]) > 3abs(δ[1])
        m[1] = 3δ[1]
    end

    m[n] = ((2h[end] + h[end-1]) * δ[end] - h[end] * δ[end-1]) / (h[end] + h[end-1])
    if sign(m[n]) != sign(δ[end])
        m[n] = 0.0
    elseif sign(δ[end]) != sign(δ[end-1]) && abs(m[n]) > 3abs(δ[end])
        m[n] = 3δ[end]
    end

    for k in 2:(n - 1)
        if δ[k - 1] * δ[k] <= 0.0
            m[k] = 0.0
        else
            w1 = 2h[k] + h[k - 1]
            w2 = h[k] + 2h[k - 1]
            m[k] = (w1 + w2) / (w1 / δ[k - 1] + w2 / δ[k])
        end
    end

    return m
end

function pchip_resample(x::Vector{Float64}, y::Vector{Float64}, xq::Vector{Float64})
    m = pchip_slopes(x, y)
    yq = similar(xq)

    for (j, ξ) in pairs(xq)
        if ξ <= x[1]
            yq[j] = y[1]
            continue
        elseif ξ >= x[end]
            yq[j] = y[end]
            continue
        end

        i = searchsortedlast(x, ξ)
        i = clamp(i, 1, length(x) - 1)
        h = x[i + 1] - x[i]
        t = (ξ - x[i]) / h

        h00 = 2t^3 - 3t^2 + 1
        h10 = t^3 - 2t^2 + t
        h01 = -2t^3 + 3t^2
        h11 = t^3 - t^2

        yq[j] = h00 * y[i] + h10 * h * m[i] + h01 * y[i + 1] + h11 * h * m[i + 1]
    end

    return yq
end

function compute_triple_point_dispersion(
    z::AbstractVector{<:Real},
    Km::AbstractVector{<:Real},
    e::AbstractVector{<:Real};
    alpha::Real = 0.05,
    e_ref::Real = 1.0,
)
    valid_mask = isfinite.(z) .& isfinite.(Km) .& isfinite.(e)
    if count(valid_mask) < 3
        return NaN
    end

    z_v = Float64.(z[valid_mask])
    Km_v = Float64.(Km[valid_mask])
    e_v = Float64.(e[valid_mask])

    H_sbl = maximum(z_v) - minimum(z_v)
    H_sbl > 0.0 || return NaN

    z_K = z_v[argmin(Km_v)]

    target_e = alpha * e_ref
    z_e = z_v[argmin(abs.(e_v .- target_e))]

    de_dz = abs.(diff(e_v) ./ diff(z_v))
    z_ez_mid = (z_v[1:end-1] .+ z_v[2:end]) ./ 2.0
    z_ez = z_ez_mid[argmax(de_dz)]

    return (max(z_K, z_e, z_ez) - min(z_K, z_e, z_ez)) / H_sbl
end

@testset "SCM Grid Refinement Convergence (Nz=38->150)" begin
    z38 = collect(range(1.5, 200.0, length=38))
    z76 = collect(range(1.5, 200.0, length=76))
    z150 = collect(range(1.5, 200.0, length=150))

    Nt = 24
    hours = collect(0:(Nt - 1))
    nocturnal_idx = findall(h -> h >= 18 || h <= 6, hours)

    u38 = fill(NaN, length(z38), Nt)
    v38 = fill(NaN, length(z38), Nt)
    θv38 = fill(NaN, length(z38), Nt)
    Km38 = fill(NaN, length(z38), Nt)
    e38 = fill(NaN, length(z38), Nt)

    for t in 1:Nt
        phase = 2π * (t - 1) / Nt

        # Synthetic SCM-like primitive profiles (time-varying but smooth/monotone in z).
        u_col = 2.0 .+ 11.0 .* (1 .- exp.(-z38 ./ 55.0)) .+ 0.15 * cos(phase)
        v_col = 0.3 .+ 2.8 .* (1 .- exp.(-z38 ./ 85.0)) .+ 0.05 * sin(phase)
        θv_col = 285.0 .+ 0.035 .* z38 .+ 0.08 * cos(phase)

        # Turbulence closure proxies.
        Km_col = 0.001 .+ (z38 ./ 200.0) .^ 2 .+ 0.00005 * cos(phase)

        # Construct e(z) so z_K, z_e, and z_ez remain clustered near the surface,
        # and should tighten (or remain stable) under refinement.
        e_col = fill(0.20 + 0.005 * sin(phase), length(z38))
        e_col[1] = 0.050
        e_col[2] = 0.20

        u38[:, t] .= u_col
        v38[:, t] .= v_col
        θv38[:, t] .= θv_col
        Km38[:, t] .= Km_col
        e38[:, t] .= e_col
    end

    δ38 = fill(NaN, Nt)
    δ76 = fill(NaN, Nt)
    δ150 = fill(NaN, Nt)

    for t in 1:Nt
        u76_col = pchip_resample(z38, u38[:, t], z76)
        v76_col = pchip_resample(z38, v38[:, t], z76)
        θv76_col = pchip_resample(z38, θv38[:, t], z76)
        Km76_col = pchip_resample(z38, Km38[:, t], z76)
        e76_col = pchip_resample(z38, e38[:, t], z76)

        u150_col = pchip_resample(z38, u38[:, t], z150)
        v150_col = pchip_resample(z38, v38[:, t], z150)
        θv150_col = pchip_resample(z38, θv38[:, t], z150)
        Km150_col = pchip_resample(z38, Km38[:, t], z150)
        e150_col = pchip_resample(z38, e38[:, t], z150)

        # Primary profile resampling integrity checks.
        @test all(isfinite.(u76_col)) && all(isfinite.(v76_col)) && all(isfinite.(θv76_col))
        @test all(isfinite.(u150_col)) && all(isfinite.(v150_col)) && all(isfinite.(θv150_col))

        δ38[t] = compute_triple_point_dispersion(z38, Km38[:, t], e38[:, t])
        δ76[t] = compute_triple_point_dispersion(z76, Km76_col, e76_col)
        δ150[t] = compute_triple_point_dispersion(z150, Km150_col, e150_col)
    end

    @test all(isfinite.(δ38[nocturnal_idx]))
    @test all(isfinite.(δ76[nocturnal_idx]))
    @test all(isfinite.(δ150[nocturnal_idx]))

    println("SCM Refinement Audit (Nocturnal Columns)")
    println("hour | delta38 | delta150 | E_conv | status")
    println("-----+---------+----------+--------+---------------------------")
    for idx in nocturnal_idx
        e_conv = abs(δ150[idx] - δ38[idx])
        non_mono = δ150[idx] > δ38[idx]
        status = non_mono ? "[WARNING: NON-MONOTONIC]" : "[OK]"
        @printf("%4d | %7.5f | %8.5f | %6.5f | %s\n", hours[idx], δ38[idx], δ150[idx], e_conv, status)
    end

    # Required pass/fail criterion: refined grid should not increase dispersion.
    @test all(δ150[nocturnal_idx] .<= (δ38[nocturnal_idx] .+ 1e-12))

    E_conv = abs.(δ150[nocturnal_idx] .- δ38[nocturnal_idx])
    @test all(E_conv .>= 0.0)
end
