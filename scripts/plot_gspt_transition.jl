using Pkg
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(PROJECT_ROOT)
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
using SBLToolkit
using Printf, Plots, Dates

const GSPT_DEBUG_AUDIT = get(ENV, "GSPT_DEBUG_AUDIT", "1") == "1"
const GSPT_NOCTURNAL_ONLY = get(ENV, "GSPT_NOCTURNAL_ONLY", "1") == "1"
const GSPT_WRITE_RAW_DIAGNOSTIC = get(ENV, "GSPT_WRITE_RAW_DIAGNOSTIC", "1") == "1"
const GSPT_SYNTHESIZE_MISSING_FLUXES = get(ENV, "GSPT_SYNTHESIZE_MISSING_FLUXES", "0") == "1"
const GSPT_HSBL_DIAGNOSTIC = parse(Float64, get(ENV, "GSPT_HSBL_DIAGNOSTIC", "200.0"))

# Override campaign sign convention when source data encodes cooling as positive flux.
const COOLING_FLUX_SIGN_BY_CAMPAIGN = Dict(
    "cases99" => :negative,
    "gabls3" => :negative,
    "floss" => :negative,
    "bllast" => :negative,
    "sheba" => :negative,
)

function run_campaign_gspt_transition(nc_path::String, campaign_name::String)
    output_dir = joinpath(@__DIR__, "..", "reports", "generated", "gspt_phase2")
    mkpath(output_dir)
    save_path = joinpath(output_dir, "$(campaign_name)_gspt_transition_heatmap.png")
    raw_save_path = joinpath(output_dir, "$(campaign_name)_gspt_transition_heatmap_raw.png")

    cooling_flux_sign = get(COOLING_FLUX_SIGN_BY_CAMPAIGN, campaign_name, :negative)

    println("Processing $(campaign_name) time-series from: $(nc_path)")

    results = ingest_netcdf_gspt(
        nc_path;
        S2_min=1e-3,
        σ_u=0.05,
        σ_v=0.05,
        σ_th=0.05,
        cooling_flux_sign=cooling_flux_sign,
        nocturnal_only=GSPT_NOCTURNAL_ONLY,
        debug_audit=GSPT_DEBUG_AUDIT,
        mask_ill_conditioned_in_solver=true,
        synthesize_missing_fluxes=GSPT_SYNTHESIZE_MISSING_FLUXES,
        h_sbl_diagnostic=GSPT_HSBL_DIAGNOSTIC
    )
    p = plot_gspt_transition(
        results;
        save_path=save_path,
        mask_ill_conditioned=true,
        mask_missing=true,
        plot_title="Dynamic Nocturnal GSPT Transition Surface (R_coord)"
    )

    finite_primary = count(isfinite, results.R_coord)
    if GSPT_WRITE_RAW_DIAGNOSTIC && finite_primary == 0
        @warn "Primary GSPT output is fully masked/NaN for $(campaign_name). Running raw fallback with nocturnal_only=false and no masks."
        raw_results = ingest_netcdf_gspt(
            nc_path;
            S2_min=1e-3,
            σ_u=0.05,
            σ_v=0.05,
            σ_th=0.05,
            cooling_flux_sign=cooling_flux_sign,
            nocturnal_only=false,
            debug_audit=true,
            mask_ill_conditioned_in_solver=false,
            synthesize_missing_fluxes=true,
            h_sbl_diagnostic=GSPT_HSBL_DIAGNOSTIC
        )
        plot_gspt_transition(
            raw_results;
            save_path=raw_save_path,
            mask_ill_conditioned=false,
            mask_missing=false,
            plot_title="Dynamic GSPT Transition Surface (R_coord) [Raw, All Timesteps]"
        )
    end

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
        nc_path = nc_files[1]
        campaign_name = campaign
        if ispath(nc_path)
            run_campaign_gspt_transition(nc_path, campaign_name)
        else
            @info "Skipping unmounted campaign dataset: $campaign_name"
        end
    else
        @warn "No .nc files found in campaign directory: $(camp_dir)"
    end
end