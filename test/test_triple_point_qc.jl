using Test
using Statistics

"""
    compute_triple_point_dispersion(z, Km, e; alpha=0.05, e_ref=1.0)

Computes non-dimensional triple-point dispersion (δ_TP) for a single vertical profile column:
1. Diffusivity cutoff height z_K = argmin(K_m)
2. TKE level-set height z_e where e(z_e) ≈ α * e_ref
3. TKE gradient extremum height z_e_z = argmax(|∂e/∂z|)
"""
function compute_triple_point_dispersion(
    z::AbstractVector{<:Real},
    Km::AbstractVector{<:Real},
    e::AbstractVector{<:Real};
    alpha::Real = 0.05,
    e_ref::Real = 1.0
)
    # Filter out missing/NaN values prior to computation
    valid_mask = isfinite.(Km) .& isfinite.(e) .& isfinite.(z)
    if count(valid_mask) < 3
        return NaN
    end

    z_v  = z[valid_mask]
    Km_v = Km[valid_mask]
    e_v  = e[valid_mask]

    # 1. Diffusivity cutoff height
    z_K = z_v[argmin(Km_v)]

    # 2. TKE level-set height
    target_e = alpha * e_ref
    z_e = z_v[argmin(abs.(e_v .- target_e))]

    # 3. TKE gradient extremum height
    de_dz = abs.(diff(e_v) ./ diff(z_v))
    z_ez_mid = (z_v[1:end-1] .+ z_v[2:end]) ./ 2.0
    z_ez = z_ez_mid[argmax(de_dz)]

    # Domain height scaling factor
    H_sbl = maximum(z_v) - minimum(z_v)
    H_sbl > 0 || return NaN

    # Non-Dimensional Triple-Point Dispersion
    return (max(z_K, z_e, z_ez) - min(z_K, z_e, z_ez)) / H_sbl
end

"""
    filter_ingested_matrix_qc(z, Km_mat, e_mat; max_dispersion=0.10)

Applies δ_TP quality-control gating to an ingested profile matrix (N_z × N_t).
Returns a boolean vector indicating valid physical columns.
"""
function filter_ingested_matrix_qc(
    z::AbstractVector{<:Real},
    Km_mat::AbstractMatrix{<:Real},
    e_mat::AbstractMatrix{<:Real};
    max_dispersion::Float64 = 0.10
)
    Nt = size(Km_mat, 2)
    qc_mask = fill(false, Nt)

    for t in 1:Nt
        δ_TP = compute_triple_point_dispersion(z, view(Km_mat, :, t), view(e_mat, :, t))
        qc_mask[t] = isfinite(δ_TP) && (δ_TP < max_dispersion)
    end

    return qc_mask
end

@testset "GSPT Ingestion Quality-Control Gate (δ_TP < 0.10)" begin
    # Standard GABLS3/CASES-99 grid vertical resolution
    z = collect(range(1.5, 200.0, length=38))
    Nt = 20

    @testset "Single Profile Metric Tests" begin
        # Synthetic aligned nocturnal collapse profile
        Km_phys = 1.0 .+ (z ./ 200.0).^2
        Km_phys[1:4] .= 0.001  # Minimum near surface

        # Keep all three triple-point markers near the surface so δ_TP is small.
        # target_e = alpha * e_ref = 0.05 (default), and a sharp early jump forces
        # |de/dz| maximum to occur near the bottom levels.
        e_phys = fill(0.20, length(z))
        e_phys[1] = 0.05
        e_phys[2] = 0.051

        δ_TP = compute_triple_point_dispersion(z, Km_phys, e_phys)

        @test isfinite(δ_TP)
        @test δ_TP < 0.10
    end

    @testset "Ingested Matrix Column Gating" begin
        # Construct synthetic profile matrix: 15 physical columns, 5 noisy/divergent columns
        Km_good = 0.001 .+ (z ./ 200.0).^2
        e_good = fill(0.20, length(z))
        e_good[1] = 0.05
        e_good[2] = 0.051

        # Deliberately divergent profiles force large triple-point spread.
        Km_bad = 0.001 .+ ((maximum(z) .- z) ./ 200.0).^2
        e_bad = 0.05 .+ 0.20 .* exp.(-(z .- 180.0).^2 ./ (2.0 * 6.0^2))

        Km_mat = hcat([Km_good for _ in 1:15]..., [Km_bad for _ in 1:5]...)
        e_mat  = hcat([e_good for _ in 1:15]..., [e_bad for _ in 1:5]...)

        qc_mask = filter_ingested_matrix_qc(z, Km_mat, e_mat; max_dispersion=0.10)

        # Confirm exact column rejection behavior
        @test length(qc_mask) == Nt
        @test all(qc_mask[1:15])       # Aligned profiles pass
        @test !any(qc_mask[16:end])    # Unphysical noisy profiles dropped
    end
end