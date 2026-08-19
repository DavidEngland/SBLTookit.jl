using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using SBLToolkit
using JLD2

if length(ARGS) < 3
    println("Usage: julia process_campaign.jl <site_name> <input_nc> <output_dir>")
    exit(1)
end

site_name = ARGS[1]
input_nc  = ARGS[2]
output_dir = ARGS[3]

mkpath(output_dir)

ds = ingest_nc(input_nc, site_name)
U_clean = clean_gaps(ds.U)

V, PC, Ubar = modal_decomposition(U_clean)

dt_seconds = Float64(Dates.value(ds.time[2] - ds.time[1])) / 1000.0
PC2_LLJ, PC2_IGW = separate_jet_and_wave(PC[2, :], dt_seconds)

U_llj, z_llj, u_max = track_jet_core(ds.heights, V[:, 2], PC2_LLJ, Ubar)

out_filename = joinpath(output_dir, "$(site_name)_processed.jld2")
jldsave(out_filename;
    site=site_name,
    time=ds.time,
    heights=ds.heights,
    V=V,
    PC=PC,
    PC2_LLJ=PC2_LLJ,
    PC2_IGW=PC2_IGW,
    z_llj=z_llj,
    u_max=u_max,
    qc_mask=ds.qc_mask
)

println("Successfully processed and saved: $out_filename")