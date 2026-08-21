#!/usr/bin/env julia
# scripts/audit_seb_nc.jl
#
# Scan NetCDF files for variables potentially relevant to the
# surface energy budget (SEB), using variable names and metadata only.

using NCDatasets

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

const NC_EXTENSIONS = Set([".nc", ".nc4", ".cdf"])

const SEB_PATTERNS = [
    # Ground / surface state
    "soil",
    "snow",
    "ice",
    "thermistor",
    "ground.?heat",
    "skin",
    "surf_temp",
    "surface_temperature",

    # Surface flux terminology
    "surface_heat_flux",
    "sensible_heat_flux",
    "latent_heat_flux",
    "surface_energy",
    "energy_flux",

    # Radiation terminology
    "rad_flux",
    "net_rad",
    "netrad",
    "downward_shortwave",
    "upward_shortwave",
    "downward_longwave",
    "upward_longwave",

    # Common model / observational abbreviations
    "^g$",
    "^g_",
    "^ts$",
    "^ts_",
    "^swc",
    "^shf$",
    "^lhf$",
    "^l hf$",
    "tskin",
    "tsfc",
    "t_sfc",
    "tice",
    "tsnw",
    "tsnow",
    "qsub",
    "^rn$",
    "^lwd$",
    "^lwu$",
    "^swd$",
    "^swu$",
]

const SEB_REGEX = Regex(join(SEB_PATTERNS, "|"), "i")

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

@inline matches_seb(name::AbstractString) =
    occursin(SEB_REGEX, name)

function default_root()
    sibling = normpath(
        joinpath(
            @__DIR__,
            "..",
            "..",
            "SpectralBL-Analytics",
            "data",
        ),
    )

    if isdir(sibling)
        return sibling
    end

    # Explicitly report the fallback rather than silently hiding it.
    @warn "Expected SpectralBL-Analytics/data directory not found; using current working directory" \
          expected_root=sibling fallback=pwd()

    return pwd()
end

function find_nc_files(root::AbstractString)
    files = String[]

    for (dirpath, _, filenames) in walkdir(root)
        for fn in filenames
            ext = lowercase(splitext(fn)[2])

            if ext in NC_EXTENSIONS
                push!(files, joinpath(dirpath, fn))
            end
        end
    end

    return sort!(files)
end

function scan_file(path::AbstractString)
    hits = Tuple{String,String,String,String}[]
    # (variable_name, long_name, standard_name, units)

    try
        NCDataset(path, "r") do ds
            for (vname, var) in ds
                vstr = String(vname)

                long_name =
                    string(get(var.attrib, "long_name", ""))

                standard_name =
                    string(get(var.attrib, "standard_name", ""))

                units =
                    string(get(var.attrib, "units", ""))

                # Search variable name and metadata independently.
                if matches_seb(vstr) ||
                   matches_seb(long_name) ||
                   matches_seb(standard_name)

                    push!(
                        hits,
                        (vstr, long_name, standard_name, units),
                    )
                end
            end
        end

    catch e
        @warn "Could not read NetCDF file" \
              path \
              exception=(e, catch_backtrace())
    end

    return hits
end

# ----------------------------------------------------------------------
# Main audit
# ----------------------------------------------------------------------

function main(args::Vector{String} = ARGS)

    root =
        length(args) >= 1 ?
        normpath(args[1]) :
        default_root()

    if !isdir(root)
        println(stderr, "Root directory not found: $root")
        return 1
    end

    println("Scanning for NetCDF files under:")
    println("  $root")
    println()

    nc_files = find_nc_files(root)

    println("Found $(length(nc_files)) NetCDF file(s).")
    println()

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

            label =
                isempty(long_name) ?
                vname :
                "$vname ($long_name)"

            if !isempty(standard_name)
                label = "$label {standard_name=$standard_name}"
            end

            if !isempty(units)
                label = "$label [$units]"
            end

            println("    ", label)
        end

        println()
    end

    if !any_hits
        println("  (no SEB-like variable names or metadata found)")
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end