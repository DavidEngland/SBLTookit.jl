#!/usr/bin/env julia
# scripts/audit_seb_nc.jl — Scan NetCDF files for surface energy budget (SEB) fields.

using NCDatasets



# Single compiled regex for faster attribute scanning
const SEB_PATTERNS = [
    "soil", "snow", "ice", "thermistor", "ground.?heat", "skin", "surface",
    "radiation", "moist", "^g$", "^g_", "^ts$", "^ts_", "^swc", "^shf$",
    "tskin", "tsfc", "t_sfc", "tice", "tsnw", "tsnow", "qsub", "^rn$",
    "netrad", "^lwd$", "^lwu$", "^swd$", "^swu$"
]
const SEB_REGEX = Regex(join(SEB_PATTERNS, "|"), "i")

matches_seb(name::AbstractString) = occursin(SEB_REGEX, name)

function default_root()
    sibling = normpath(joinpath(@__DIR__, "..", "..", "SpectralBL-Analytics", "data"))
    return isdir(sibling) ? sibling : pwd()
end

function find_nc_files(root::AbstractString)
    files = String[]
    for (dirpath, _, filenames) in walkdir(root)
        for fn in filenames
            if lowercase(splitext(fn)[2]) == ".nc"
                push!(files, joinpath(dirpath, fn))
            end
        end
    end
    return sort!(files)
end

function scan_file(path::AbstractString)
    hits = Tuple{String,String,String,String}[] # (vname, long_name, standard_name, units)
    try
        NCDataset(path, "r") do ds
            for (vname, var) in ds
                vstr = String(vname)
                long_name = string(get(var.attrib, "long_name", ""))
                standard_name = string(get(var.attrib, "standard_name", ""))
                units = string(get(var.attrib, "units", ""))

                if any(matches_seb, (vstr, long_name, standard_name))
                    push!(hits, (vstr, long_name, standard_name, units))
                end
            end
        end
    catch e
        @warn "Could not read NetCDF file" path exception = e
    end
    return hits
end

function main(args::Vector{String} = ARGS)
    root = length(args) >= 1 ? normpath(args[1]) : default_root()
    if !isdir(root)
        println(stderr, "Root directory not found: $root")
        return 1
    end

    println("Scanning for .nc files under: $root")
    nc_files = find_nc_files(root)
    println("Found $(length(nc_files)) NetCDF file(s).\n")

    println("============================================")
    println(" NetCDF variables matching SEB terms")
    println("============================================")

    any_hits = false
    for f in nc_files
        hits = scan_file(f)
        isempty(hits) && continue
        any_hits = true

        display_path = relpath(f, root)
        println("-- ", display_path, " --")
        for (vname, long_name, standard_name, units) in hits
            label = isempty(long_name) ? vname : "$vname ($long_name)"
            label = isempty(standard_name) ? label : "$label {standard_name=$standard_name}"
            label = isempty(units) ? label : "$label [$units]"
            println("    ", label)
        end
        println()
    end

    any_hits || println("  (no SEB-like variable names found)")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end