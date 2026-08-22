using Dates
using Glob

struct CampaignFileSpec
    path::String
    date::Date
    adapter::Symbol
end

function discover_campaign_files(config::Dict)::Vector{CampaignFileSpec}
    root = config["root_dir"]
    pattern = Regex(config["file_pattern"])
    recursive = get(config, "recursive", false)
    adapter = Symbol(config["adapter"])

    # Fast directory walk
    files = String[]
    if recursive
        for (dirpath, _, filenames) in walkdir(root)
            append!(files, [joinpath(dirpath, f) for f in filenames])
        end
    else
        files = readdir(root, join=true)
    end

    matched_files = CampaignFileSpec[]

    for file_path in files
        m = match(pattern, basename(file_path))
        isnothing(m) && continue

        # Extract date components from named regex capture groups
        year = if haskey(m, :year)
            parse(Int, m[:year])
        elseif haskey(m, :y2)
            parse(Int, m[:y2]) + get(config, "year_offset", 2000)
        else
            get(config, "base_year", 2026)
        end

        month = parse(Int, m[:month])
        day = parse(Int, m[:day])

        push!(matched_files, CampaignFileSpec(file_path, Date(year, month, day), adapter))
    end

    # Sort files chronologically
    sort!(matched_files, by=x -> x.date)
    return matched_files
end