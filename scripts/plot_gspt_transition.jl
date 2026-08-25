using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
using SBLToolkit
using Printf, Plots, Dates

function run_campaign_gspt_transition(nc_path::String, campaign_name::String)
    output_dir = joinpath(@__DIR__, "..", "reports", "generated", "gspt_phase2")
    mkpath(output_dir)
    save_path = joinpath(output_dir, "$(campaign_name)_gspt_transition_heatmap.png")

    println("Processing $(campaign_name) time-series from: $(nc_path)")

    results = ingest_netcdf_gspt(nc_path; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
    p = plot_gspt_transition(results; save_path=save_path)

    println("Saved $(campaign_name) transition heatmap to $(save_path)\n")
    return results
end

raw_dir = joinpath(@__DIR__, "..", "data", "raw")
campaigns = ["cases99", "gabls3", "floss", "bllast", "sheba"]

for campaign in campaigns
    camp_dir = joinpath(raw_dir, campaign)

    if !ispath(camp_dir)
        @warn "Campaign directory path or symlink target does not exist: $(camp_dir)"
        continue
    end

    nc_files = String[]
    for (root, _, files) in walkdir(camp_dir; follow_symlinks=true)
        for file in files
            if endswith(lowercase(file), ".nc")
                push!(nc_files, joinpath(root, file))
            end
        end
    end

    if !isempty(nc_files)
        # Prioritize root-level or uniform processing files over raw nested subfolders
        sort!(nc_files, by=path -> contains(path, "/raw/") ? 2 : 1)
        run_campaign_gspt_transition(nc_files[1], campaign)
    else
        @warn "No .nc files found in campaign directory: $(camp_dir)"
    end
end