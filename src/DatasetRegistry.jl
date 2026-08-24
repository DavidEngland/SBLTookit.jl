module DatasetRegistry

using JSON
using DataFrames

export discover_campaign_trajectories, resolve_dataset_path

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const ANALYTICS_ROOT = normpath(joinpath(PROJECT_ROOT, "..", "SpectralBL-Analytics"))

function resolve_manifest_path()
    local_path = joinpath(PROJECT_ROOT, "data", "datasets.json")
    isfile(local_path) && return local_path

    analytics_path = joinpath(ANALYTICS_ROOT, "data", "datasets.json")
    isfile(analytics_path) && return analytics_path

    return local_path
end

function normalize_campaign_name(raw::String)
    clean = replace(lowercase(raw), "_" => "")
    if clean in ("cases99", "cases")
        return "CASES-99"
    elseif startswith(clean, "floss")
        return uppercase(raw)
    elseif clean == "bllast"
        return "BLLAST"
    elseif clean == "sheba"
        return "SHEBA"
    elseif clean in ("gabls3", "gabs3")
        return "GABLS3"
    end
    return uppercase(raw)
end

function fallback_trajectory_scan()
    candidate_dirs = [
        joinpath(PROJECT_ROOT, "data", "drafts", "trajectories"),
        joinpath(ANALYTICS_ROOT, "data", "drafts", "trajectories")
    ]
    map = Dict{String,String}()
    for drafts_dir in candidate_dirs
        !isdir(drafts_dir) && continue
        for f in readdir(drafts_dir; join=true)
            if endswith(f, ".csv") && startswith(basename(f), "trajectory_")
                raw = replace(basename(f), r"^trajectory_|\.csv$" => "")
                map[normalize_campaign_name(raw)] = f
            end
        end
        !isempty(map) && break
    end
    return map
end

"""
    discover_campaign_trajectories([manifest_path::String]) -> Dict{String, String}

Parses datasets.json (supporting both Array catalog and Dict schema formats)
and returns a map of campaign display names to resolved trajectory CSV paths.
"""
function discover_campaign_trajectories(manifest_path::String)
    if !isfile(manifest_path)
        return fallback_trajectory_scan()
    end

    manifest = JSON.parsefile(manifest_path)
    datasets = get(manifest, "datasets", Dict())
    campaign_map = Dict{String,String}()

    # Catalog Array Schema
    if datasets isa Vector
        for item in datasets
            status = get(item, "status", "")
            data_assets = get(item, "data_assets", Dict())
            traj_rel = get(data_assets, "trajectory_csv", "")

            if (status == "Ingested" || !isempty(traj_rel)) && !isempty(traj_rel)
                p_local = joinpath(PROJECT_ROOT, traj_rel)
                p_analytics = joinpath(ANALYTICS_ROOT, traj_rel)

                real_path = isfile(p_local) ? p_local : (isfile(p_analytics) ? p_analytics : nothing)
                if real_path !== nothing
                    id_raw = get(item, "id", get(item, "name", ""))
                    display_name = normalize_campaign_name(id_raw)
                    campaign_map[display_name] = real_path
                end
            end
        end

        # Legacy Dict Schema
    elseif datasets isa Dict
        trajectory_config = get(datasets, "pipeline_draft_trajectories", nothing)
        if trajectory_config !== nothing && get(trajectory_config, "ingest", false)
            base_rel = get(trajectory_config, "base_dir", "")
            base_dir = isdir(joinpath(PROJECT_ROOT, base_rel)) ?
                       joinpath(PROJECT_ROOT, base_rel) :
                       joinpath(ANALYTICS_ROOT, base_rel)

            pattern_str = get(trajectory_config, "file_pattern", "")
            if !isempty(pattern_str)
                regex = Regex(pattern_str)
                if isdir(base_dir)
                    for file in readdir(base_dir; join=false)
                        m = match(regex, file)
                        if m !== nothing && haskey(m, :campaign)
                            raw_campaign = m[:campaign]
                            display_name = normalize_campaign_name(raw_campaign)
                            campaign_map[display_name] = joinpath(base_dir, file)
                        end
                    end
                end
            end
        end
    end

    return isempty(campaign_map) ? fallback_trajectory_scan() : campaign_map
end

discover_campaign_trajectories() = discover_campaign_trajectories(resolve_manifest_path())

function resolve_dataset_path(dataset_key::String, manifest_path::String)
    if !isfile(manifest_path)
        return nothing
    end
    manifest = JSON.parsefile(manifest_path)
    datasets = get(manifest, "datasets", Dict())

    if datasets isa Dict
        ds = get(datasets, dataset_key, nothing)
        ds === nothing && return nothing
        return joinpath(PROJECT_ROOT, get(ds, "base_dir", ""))
    elseif datasets isa Vector
        for item in datasets
            if get(item, "id", "") == dataset_key
                data_assets = get(item, "data_assets", Dict())
                traj = get(data_assets, "trajectory_csv", "")
                if !isempty(traj)
                    p_local = joinpath(PROJECT_ROOT, traj)
                    p_analytics = joinpath(ANALYTICS_ROOT, traj)
                    return isfile(p_local) ? dirname(p_local) : dirname(p_analytics)
                end
            end
        end
    end
    return nothing
end

resolve_dataset_path(dataset_key::String) = resolve_dataset_path(dataset_key, resolve_manifest_path())

end # module DatasetRegistry