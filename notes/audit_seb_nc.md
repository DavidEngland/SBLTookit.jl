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




—-

Broad Pattern Matches (High False-Positive Rate)

Patterns like "surface", "radiation", and "moist" will trigger on non-SEB atmospheric and coordinate fields (e.g., surface_pressure, surface_altitude, radiation_scheme_id, or atmospheric specific_humidity). If you want to isolate energy fluxes and ground state variables, narrow these tokens:

Julia

# Replace broad tokens: "surface", "radiation", "moist"
# With specific SEB terms:
const SEB_PATTERNS = [
    "soil", "snow", "ice", "thermistor", "ground.?heat", "skin",
    "surf_temp", "surface_temperature", "surface_heat_flux",
    "rad_flux", "net_rad", "downward_shortwave", "upward_longwave",
    "^g$", "^g_", "^ts$", "^ts_", "^swc", "^shf$",
    "tskin", "tsfc", "t_sfc", "tice", "tsnw", "tsnow", "qsub", "^rn$",
    "netrad", "^lwd$", "^lwu$", "^swd$", "^swu$"
]
Missing File Extensions

splitext(fn)[2] currently only checks for .nc. NetCDF files frequently use .nc4 or .cdf extensions.

Julia

# Expand extension check
const NC_EXTENSIONS = Set([".nc", ".nc4", ".cdf"])
if lowercase(splitext(fn)[2]) in NC_EXTENSIONS
Hardcoded Relative Path Dependency

default_root() relies on a rigid sibling directory structure (../../SpectralBL-Analytics/data). If running outside that specific workspace layout, it defaults silently to pwd().

Key Strengths

* Metadata-Only Processing: Inspecting ds keys and var.attrib without loading tensor data keeps memory footprint minimal and runtime fast across large datasets.
* Global Regex Pre-compilation: Top-level instantiation of SEB_REGEX avoids redundant string concatenation and compile overhead per file.
* Graceful Exception Handling: The try...catch block around NCDataset prevents corrupted files from crashing the audit sweep.

—-

Yes. The two proposed corrections are worthwhile, and I would make one additional change: separate variable-name matching from attribute matching. Right now the same broad regex is applied to vname, long_name, and standard_name, so even a relatively narrow pattern can produce false positives through metadata.
Recommended changes
1. Use .nc, .nc4, and .cdf
    * Definitely appropriate for a dataset audit.
    * A Set makes the intent explicit.
2. Remove broad surface, radiation, and moist
    * Agreed. Those are much too permissive for an SEB audit.
    * In particular, "moist" can hit atmospheric moisture variables that have nothing to do with surface energy balance.
3. Use SEB-specific terminology
    * Keep soil, snow, ice, skin, thermistor, etc.
    * Add explicit radiation/flux terms.
    * Keep conventional abbreviations such as shf, g, rn, lwd, etc.
4. Improve the root fallback
    * Silently falling back to pwd() can cause an apparently successful audit of the wrong directory.
    * Better behavior is to search the expected location, then explicitly report that the working directory is being used.
5. One important technical issue: matches_seb is currently called with a tuple:any(matches_seb, (vstr, long_name, standard_name))This is valid because any(f, itr) applies f elementwise, but the name matches_seb actually expects a string. That’s fine, just slightly opaque. I would make the matching logic explicit.
I would use this version
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
There is one typo I would correct before using that verbatim: "^l hf$" should not be in the pattern list. It appears to be an accidental space-containing variant and has no useful purpose. Remove it.
More importantly, for your SHEBA/SEB audit, I would consider taking the next step beyond simple regex matching: classify hits into surface temperature, radiative flux, turbulent flux, ground/soil heat, snow/ice state, and hydrological/subli­mation terms. That would turn this from a file-name scanner into a genuinely useful SEB field inventory for determining whether the raw datasets contain enough information to diagnose the SCM surface-energy-budget instability.
