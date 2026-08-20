#!/usr/bin/env julia
# scripts/audit_campaign_data.jl — Audit centralized campaign data across NetCDF, CSV/ASCII, and documentation.
# Audit centralized campaign data across NetCDF, CSV/ASCII, and documentation.
#
# Usage:
#   julia --project=. scripts/audit_campaign_data.jl [root_dir]
#
# If root_dir is omitted, GeoABL.get_data_dir() resolves GEOABL_DATA_DIR,
# SpectralBL-Analytics/data, or the local data fallback.

using CSV
using GeoABL: get_data_dir
using NCDatasets

const SEB_PATTERN_STRINGS = [
    "soil", "snow", "ice", "thermistor", "ground.?heat", "skin", "surface",
    "radiation", "moist", "^g$", "^g_", "^ts$", "^ts_", "^swc", "^shf$",
    "tskin", "tsfc", "t_sfc", "tice", "tsnw", "tsnow", "qsub", "^rn$",
    "netrad", "^lwd$", "^lwu$", "^swd$", "^swu$",
]
const SEB_REGEX = Regex(join(SEB_PATTERN_STRINGS, "|"), "i")

const DOCUMENT_EXTENSIONS = Set([".md", ".pdf", ".rst"])
const TABULAR_EXTENSIONS = Set([".csv", ".dat", ".asc", ".tsv", ".txt"])
const DOCUMENT_NAME_PATTERN =
    r"readme|metadata|schema|source|description|overview|background"i

matches_seb(value::AbstractString) = occursin(SEB_REGEX, value)

function classify_files(root_dir::AbstractString)
    netcdf_files = String[]
    tabular_files = String[]
    documentation_files = String[]

    for (dirpath, _, filenames) in walkdir(root_dir)
        for filename in filenames
            path = joinpath(dirpath, filename)
            extension = lowercase(splitext(filename)[2])

            if extension == ".nc"
                push!(netcdf_files, path)
            elseif extension in DOCUMENT_EXTENSIONS ||
                   (extension == ".txt" && occursin(DOCUMENT_NAME_PATTERN, filename))
                push!(documentation_files, path)
            elseif extension in TABULAR_EXTENSIONS
                push!(tabular_files, path)
            end
        end
    end

    return sort!(netcdf_files), sort!(tabular_files), sort!(documentation_files)
end

function audit_netcdf(files::Vector{String}, root_dir::AbstractString)
    println("--- NetCDF SEB Variables ---")
    found = false

    for path in files
        try
            hits = String[]
            NCDataset(path, "r") do dataset
                for (variable_name, variable) in dataset
                    name = String(variable_name)
                    long_name = string(get(variable.attrib, "long_name", ""))
                    standard_name = string(get(variable.attrib, "standard_name", ""))
                    units = string(get(variable.attrib, "units", ""))

                    if any(matches_seb, (name, long_name, standard_name))
                        description = isempty(long_name) ? name : "$name ($long_name)"
                        description =
                            isempty(standard_name) ? description :
                            "$description {standard_name=$standard_name}"
                        description = isempty(units) ? description : "$description [$units]"
                        push!(hits, description)
                    end
                end
            end

            if !isempty(hits)
                found = true
                println("  [NC] ", relpath(path, root_dir))
                foreach(hit -> println("       - ", hit), hits)
            end
        catch error
            @warn "Could not audit NetCDF file" path exception = (error, catch_backtrace())
        end
    end

    found || println("  (no matching NetCDF variables)")
    println()
end

function csv_header_tokens(path::AbstractString)
    filesize(path) == 0 && return String[]
    rows = collect(CSV.Rows(path; header = false, limit = 2, types = String))
    tokens = String[]
    for row in rows, value in row
        push!(tokens, ismissing(value) ? "" : String(value))
    end
    return tokens
end

function ascii_header_tokens(path::AbstractString)
    filesize(path) == 0 && return String[]
    lines = String[]
    open(path, "r") do io
        while !eof(io) && length(lines) < 2
            line = strip(readline(io))
            isempty(line) || push!(lines, line)
        end
    end
    isempty(lines) && return String[]
    return reduce(vcat, (split(line, r"[,;\t\s]+") for line in lines); init = String[])
end

function audit_tabular(files::Vector{String}, root_dir::AbstractString)
    println("--- Tabular (CSV/ASCII) SEB Headers ---")
    found = false

    for path in files
        try
            extension = lowercase(splitext(path)[2])
            tokens =
                extension == ".csv" ? csv_header_tokens(path) : ascii_header_tokens(path)
            hits = unique(filter(matches_seb, tokens))
            if !isempty(hits)
                found = true
                println("  [TABLE] ", relpath(path, root_dir))
                println("          Matched headers: ", join(hits, ", "))
            end
        catch error
            @warn "Could not audit tabular file" path exception = (error, catch_backtrace())
        end
    end

    found || println("  (no matching tabular headers)")
    println()
end

function audit_documentation(files::Vector{String}, root_dir::AbstractString)
    println("--- Campaign Documentation and Metadata Files ---")
    if isempty(files)
        println("  (none found)")
    else
        for path in files
            size_kib = filesize(path) / 1024
            println(
                "  [DOC] ",
                relpath(path, root_dir),
                " (",
                round(size_kib; digits = 1),
                " KiB)",
            )
        end
    end
end

function audit_data_tree(root_dir::AbstractString)
    root = normpath(root_dir)
    isdir(root) || throw(ArgumentError("Campaign data directory not found: $root"))

    netcdf_files, tabular_files, documentation_files = classify_files(root)

    println("==================================================================")
    println(" Campaign Data Audit: $root")
    println("==================================================================")
    println(
        "Summary: Found $(length(netcdf_files)) NetCDF, ",
        "$(length(tabular_files)) tabular/ASCII, and ",
        "$(length(documentation_files)) documentation file(s).\n",
    )

    audit_netcdf(netcdf_files, root)
    audit_tabular(tabular_files, root)
    audit_documentation(documentation_files, root)
    return nothing
end

function main(args::Vector{String} = ARGS)
    root = isempty(args) ? get_data_dir() : normpath(first(args))
    if !isdir(root)
        println(
            stderr,
            "Directory not found: $root. Set GEOABL_DATA_DIR or pass a path as ARGS[1].",
        )
        return 1
    end

    audit_data_tree(root)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
