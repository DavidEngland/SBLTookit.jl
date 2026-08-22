module SBLToolkit

using Dates
using DataFrames
using NCDatasets
using DSP
using LinearAlgebra
using Statistics
using Interpolations
using JSON3

struct SBLDataset
    site::String
    time::Vector{DateTime}
    heights::Vector{Float64}
    U::Matrix{Float64}         # Shape: (n_heights, n_time)
    theta::Matrix{Float64}     # Shape: (n_heights, n_time)
    u_star::Vector{Float64}
    qc_mask::BitVector
end

# -------------------------------------------------------------------
# Include and expose submodules from src/ultra
# -------------------------------------------------------------------
if isdir(joinpath(@__DIR__, "ultra"))
    include("ultra/core_types.jl")
    using .CoreTypes

    include("ultra/forcing.jl")
    using .UltraForcing

    include("ultra/stability.jl")
    using .UltraStability

    include("ultra/spectral_engine.jl")
    using .SpectralEngine

    include("ultra/spectral_reconstruction.jl")
    using .SpectralReconstruction

    # Adapters
    include("ultra/adapters/cabauw_adapter.jl")
    using .CabauwAdapter

    include("ultra/adapters/neon_adapter.jl")
    using .NEONAdapter

    include("ultra/adapters/icos_adapter.jl")
    using .ICOSAdapter

    include("ultra/adapters/ameriflux_adapter.jl")
    using .AmeriFluxAdapter

    # Ingestion Router
    include("ultra/UnifiedBLIngestion.jl")
    using .UnifiedBLIngestion

    # Re-export public API symbols at package root
    export MeteorologicalProfile, ProfileMetadata, StandardizedBLObservation
    export classify_stability, generate_scm_forcing
    export chebyshev_fingerprint, reconstruct_from_chebyshev
    export extract_temperature_observations, extract_neon_observations, upscale_sparse_icos_observation
    export load_ameriflux_registry, get_site_metadata
    export detect_network_format, ingest_boundary_layer_data, batch_spectral_fingerprints
end

# -------------------------------------------------------------------
# Core SBLToolkit Functions
# -------------------------------------------------------------------
export SBLDataset, ingest_nc, clean_gaps, modal_decomposition,
    butterworth_lowpass, separate_jet_and_wave, track_jet_core

function ingest_nc(file_path::String, site_name::String)::SBLDataset
    return NCDataset(file_path, "r") do ds
        time = ds["time"][:]
        heights = Float64.(ds["height"][:])
        U = Float64.(ds["U"][:, :])
        theta = Float64.(ds["theta"][:, :])
        u_star = Float64.(ds["u_star"][:])

        raw_qc = haskey(ds, "qc_mask") ? BitVector(ds["qc_mask"][:]) : trues(length(time))
        flatline_mask = abs.(u_star .- 0.3) .< 1e-4
        qc_mask = raw_qc .& .!flatline_mask

        SBLDataset(site_name, time, heights, U, theta, u_star, qc_mask)
    end
end

function clean_gaps(U::Matrix{Float64})::Matrix{Float64}
    U_clean = copy(U)
    n_heights, n_time = size(U)
    for z in 1:n_heights
        row = U_clean[z, :]
        nans = isnan.(row)
        if any(nans)
            valid_idx = findall(.!nans)
            itp = linear_interpolation(valid_idx, row[valid_idx], extrapolation_bc=Line())
            U_clean[z, nans] = itp(findall(nans))
        end
    end
    return U_clean
end

function modal_decomposition(U::Matrix{Float64})
    Ubar = mean(U, dims=2)
    Uanom = U .- Ubar
    F = svd(Uanom)
    V = F.U
    PC = F.S .* F.Vt
    return V, PC, vec(Ubar)
end

function butterworth_lowpass(data::Vector{Float64}, fs::Float64, fc::Float64; order::Int=4)::Vector{Float64}
    responsetype = Lowpass(fc; fs=fs)
    designmethod = Butterworth(order)
    fil = digitalfilter(responsetype, designmethod)
    return filtfilt(fil, data)
end

function separate_jet_and_wave(PC2::Vector{Float64}, dt_seconds::Float64; T_cutoff_hours::Float64=3.0)
    fs = 1.0 / dt_seconds
    fc = 1.0 / (T_cutoff_hours * 3600.0)
    PC2_LLJ = butterworth_lowpass(PC2, fs, fc)
    PC2_IGW = PC2 .- PC2_LLJ
    return PC2_LLJ, PC2_IGW
end

function track_jet_core(heights::Vector{Float64}, V2::Vector{Float64}, PC2_LLJ::Vector{Float64}, Ubar::Vector{Float64})
    n_time = length(PC2_LLJ)
    U_llj = Ubar .+ (V2 * PC2_LLJ')
    z_llj = zeros(n_time)
    u_max = zeros(n_time)

    for t in 1:n_time
        profile = U_llj[:, t]
        idx = argmax(profile)
        if 1 < idx < length(heights)
            z1, z2, z3 = heights[idx-1], heights[idx], heights[idx+1]
            u1, u2, u3 = profile[idx-1], profile[idx], profile[idx+1]
            denom = (z1 - z2)*(z1 - z3)*(z2 - z3)
            A = (z3*(u2 - u1) + z2*(u1 - u3) + z1*(u3 - u2)) / denom
            B = (z3^2*(u1 - u2) + z2^2*(u3 - u1) + z1^2*(u2 - u3)) / denom
            z_peak = -B / (2A)
            if z1 <= z_peak <= z3
                z_llj[t] = z_peak
                u_max[t] = u2 - (B^2 / (4A))
                continue
            end
        end
        z_llj[t] = heights[idx]
        u_max[t] = profile[idx]
    end
    return U_llj, z_llj, u_max
end

end # module SBLToolkit