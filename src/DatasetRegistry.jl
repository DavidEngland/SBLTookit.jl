module DatasetRegistry

using JSON
using DataFrames

export discover_campaign_trajectories, resolve_dataset_path

const PROJECT_ROOT = normpath(joinpath(@__DIR__,
".."))
const MANIFEST_PATH = joinpath(PROJECT_ROOT,
"data",
"datasets.json")

"""
    discover_campaign_trajectories() -> Dict{String, String
}

Parses datasets.json for `pipeline_draft_trajectories` and returns a map of
campaign display names (e.g.,
"CASES-99") to absolute CSV file paths.
"""
function discover_campaign_trajectories(manifest_path: :String = MANIFEST_PATH)
    if !isfile(manifest_path)
        @warn "Manifest not found at $manifest_path. Falling back to default directory scanning."
        return fallback_trajectory_scan()
    end

    manifest = JSON.parsefile(manifest_path)
    datasets = get(manifest,
"datasets", Dict())

    trajectory_config = get(datasets,
"pipeline_draft_trajectories", nothing)
    if trajectory_config === nothing || !get(trajectory_config,
"ingest", false)
        return fallback_trajectory_scan()
    end

    base_dir = joinpath(PROJECT_ROOT, trajectory_config[
    "base_dir"
])
    pattern_str = trajectory_config[
    "file_pattern"
]
    regex = Regex(pattern_str)

    campaign_map = Dict{String, String
}()

    if isdir(base_dir)
        for file in readdir(base_dir; join=false)
            m = match(regex, file)
            if m !== nothing && haskey(m,
:campaign)
                raw_campaign = m[
    :campaign
]
                display_name = normalize_campaign_name(raw_campaign)
                campaign_map[display_name
] = joinpath(base_dir, file)
            end
        end
    end

    return isempty(campaign_map) ? fallback_trajectory_scan() : campaign_map
end

function normalize_campaign_name(raw: :String)
    clean = replace(lowercase(raw),
"_" => "")
    if clean == "cases99" || clean == "cases"
        return "CASES-99"
    elseif startswith(clean,
"floss")
        return uppercase(raw) # FLOSS or FLOSS_II
    elseif clean == "bllast"
        return "BLLAST"
    elseif clean == "sheba"
        return "SHEBA"
    elseif clean == "gabls3" || clean == "gabs3"
        return "GABLS3"
    end
    return uppercase(raw)
end

function fallback_trajectory_scan()
    drafts_dir = joinpath(PROJECT_ROOT,
"data",
"drafts",
"trajectories")
    map = Dict{String, String
}()
    !isdir(drafts_dir) && return map

    for f in readdir(drafts_dir; join=true)
        if endswith(f,
".csv") && startswith(basename(f),
"trajectory_")
            raw = replace(basename(f), r"^trajectory_|\.csv$" => "")
            map[normalize_campaign_name(raw)
] = f
        end
    end
    return map
end

end # module DatasetRegistry