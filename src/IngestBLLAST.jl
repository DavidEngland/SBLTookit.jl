#!/usr/bin/env julia
# SBLTookit.jl - Ingest BLLAST 60m Mast NetCDF data with level-specific temp offsets and stationarity filtering.
# src/IngestBLLAST.jl
"""
Ingests BLLAST 60m Mast NetCDF, applying level-specific temp offsets and stationarity filtering.
"""
function ingest_bllast(file_path::String)::SBLDataset
    return NCDataset(file_path, "r") do ds
        time = ds["time"][:]
        heights = Float64.([15.0, 30.0, 45.0, 60.0]) # Primary mast levels

        U = Float64.(ds["U"][:, :])
        theta = Float64.(ds["theta"][:, :])

        # Apply BLLAST level-specific temperature offsets (K) if uncalibrated
        temp_offsets = [0.0, -0.12, +0.05, +0.18]
        theta .+= temp_offsets

        u_star = Float64.(ds["u_star"][:])

        # Stationarity flags (Flag H / LE)
        stationarity_flag = ds["flag_stationarity"][:]
        qc_mask = BitVector(stationarity_flag .== 0) # 0 = High Quality / Stationary

        SBLDataset("BLLAST", time, heights, U, theta, u_star, qc_mask)
    end
end