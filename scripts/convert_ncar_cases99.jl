#!/usr/bin/env julia
# scripts/convert_ncar_cases99.jl — Convert NCAR cases to a standard format.

using NCDatasets
using CSV
using Dates
using Logging

const DEFAULT_Q = 0.005
const P0 = 100_000.0
const R_OVER_CP = 0.2854

const TIME_ALIASES = ("time", "Time", "t", "TIME")
const HEIGHT_ALIASES = ("height", "z", "Z", "level", "ht")
const THETA_ALIASES = ("theta", "th", "pot_temp", "potential_temperature")
const TEMP_ALIASES = ("T", "temp", "temperature", "tair", "air_temperature")
const PRESSURE_ALIASES = ("P", "pres", "pressure", "air_pressure")
const U_ALIASES = ("u", "u_wind", "U")
const V_ALIASES = ("v", "v_wind", "V")
const WS_ALIASES = ("ws", "wspd", "wind_speed")
const WD_ALIASES = ("wd", "wdir", "wind_direction")
const Q_ALIASES = ("q", "qv", "specific_humidity", "mixing_ratio")

function _find_netcdf_files(root::AbstractString)
    files = String[]
    for (dir, _, names) in walkdir(root)
        for name in names
            lower = lowercase(name)
            if endswith(lower, ".nc") || endswith(lower, ".cdf")
                push!(files, joinpath(dir, name))
            end
        end
    end
    return sort(files)
end

function _resolve_var_name(
    ds::NCDataset,
    aliases::Tuple{Vararg{String}};
    required::Bool=true,
)
    for name in aliases
        if haskey(ds, name)
            return name
        end
    end
    if required
        throw(
            ArgumentError(
                "Missing required variable. Expected one of: $(join(aliases, ", "))",
            ),
        )
    end
    return nothing
end

function _parse_time_units(units::AbstractString)
    lower = lowercase(strip(units))
    scale = if startswith(lower, "second") || startswith(lower, "sec")
        1.0
    elseif startswith(lower, "minute") || startswith(lower, "min")
        60.0
    elseif startswith(lower, "hour")
        3600.0
    elseif startswith(lower, "day")
        86400.0
    else
        nothing
    end

    idx = findfirst(" since ", lower)
    idx === nothing && return nothing

    origin_str = strip(units[(idx.stop+1):end])
    origin_str = replace(origin_str, 'T' => ' ')
    origin_str = replace(origin_str, 'Z' => ' ')
    origin_str = split(origin_str, '+')[1]
    origin_str = replace(origin_str, r"\s+[+-]\d\d:?\d\d$" => "")

    formats = (
        dateformat"yyyy-mm-dd HH:MM:SS.s",
        dateformat"yyyy-mm-dd HH:MM:SS",
        dateformat"yyyy-mm-dd HH:MM",
        dateformat"yyyy-mm-dd",
    )

    parsed = nothing
    for fmt in formats
        try
            parsed = DateTime(strip(origin_str), fmt)
            break
        catch
        end
    end
    parsed === nothing && return nothing
    scale === nothing && return nothing

    return (scale=scale, origin=datetime2unix(parsed))
end

function _read_time_seconds(ds::NCDataset, time_name::String)
    t_raw = collect(ds[time_name][:])
    isempty(t_raw) && throw(ArgumentError("Empty time axis in variable '$time_name'"))

    first_val = first(t_raw)
    if first_val isa TimeType
        return Float64[datetime2unix(DateTime(t)) for t in t_raw]
    end

    t_numeric = Float64.(t_raw)
    units = try
        ds[time_name].attrib["units"]
    catch
        nothing
    end

    if units isa AbstractString
        parsed = _parse_time_units(units)
        if parsed !== nothing
            return parsed.origin .+ parsed.scale .* t_numeric
        end
    end

    return t_numeric
end

function _extract_time_series(ds::NCDataset, varname::String, nt::Int, station_index::Int)
    data = Array(ds[varname][:])
    dims = collect(dimnames(ds[varname]))
    station_count = haskey(ds.dim, "station") ? Int(ds.dim["station"]) : 1

    to_f64(x) = Float64.(coalesce.(x, NaN))

    if ndims(data) == 0
        scalar = ismissing(data[]) ? NaN : Float64(data[])
        return fill(scalar, nt)
    elseif ndims(data) == 1
        if length(data) == nt
            return to_f64(data)
        elseif length(data) == nt * station_count && ("station" in dims) && ("time" in dims)
            station_clamped = clamp(station_index, 1, station_count)
            if findfirst(==("station"), dims) == 1
                matrix = reshape(to_f64(data), station_count, nt)
                return vec(matrix[station_clamped, :])
            else
                matrix = reshape(to_f64(data), nt, station_count)
                return vec(matrix[:, station_clamped])
            end
        end
        throw(
            DimensionMismatch(
                "$varname has length $(length(data)); expected $nt or station-time flattened",
            ),
        )
    elseif ndims(data) == 2
        time_pos = findfirst(==("time"), dims)
        station_pos = findfirst(==("station"), dims)

        if time_pos == 1 && station_pos == 2
            station_clamped = clamp(station_index, 1, size(data, 2))
            return to_f64(data[:, station_clamped])
        elseif station_pos == 1 && time_pos == 2
            station_clamped = clamp(station_index, 1, size(data, 1))
            return to_f64(data[station_clamped, :])
        end

        if time_pos == 1
            return to_f64(data[:, 1])
        elseif time_pos == 2
            return to_f64(data[1, :])
        end
    end

    throw(
        DimensionMismatch(
            "Unsupported dimensions for $varname: dims=$(dims), size=$(size(data))",
        ),
    )
end

@inline _level_token_to_height(token::AbstractString) =
    parse(Float64, replace(token, "_" => "."))

function _extract_height_suffix_vars(ds::NCDataset, prefixes::Vector{String})
    out = Dict{Float64,String}()
    for name in keys(ds)
        name_str = String(name)
        m = match(r"^([A-Za-z]+)_(\d+(?:_\d+)?)m$", name_str)
        if m === nothing
            continue
        end
        if m.captures[1] in prefixes
            out[_level_token_to_height(m.captures[2])] = name_str
        end
    end
    return out
end

function _build_rows_from_wide_schema!(
    ds::NCDataset,
    t_values::Vector{Float64},
    station_index::Int,
    time_col::Vector{Float64},
    z_col::Vector{Float64},
    u_col::Vector{Float64},
    v_col::Vector{Float64},
    theta_col::Vector{Float64},
    q_col::Vector{Float64},
)
    nt = length(t_values)

    temp_map = _extract_height_suffix_vars(ds, ["T"])
    u_map = _extract_height_suffix_vars(ds, ["U", "u"])
    v_map = _extract_height_suffix_vars(ds, ["V", "v"])
    q_map = _extract_height_suffix_vars(ds, ["q", "qv"])

    p_name = _resolve_var_name(ds, PRESSURE_ALIASES; required=false)
    p_series =
        p_name === nothing ? fill(P0, nt) :
        _extract_time_series(ds, p_name, nt, station_index)
    if maximum(p_series) < 2_000.0
        p_series .*= 100.0
    end

    common_levels =
        sort(collect(intersect(keys(temp_map), intersect(keys(u_map), keys(v_map)))))
    isempty(common_levels) && throw(
        ArgumentError(
            "No common CASES-99 levels found across T_*, U_*/u_*, and V_*/v_* variables",
        ),
    )

    for z_val in common_levels
        t_name = temp_map[z_val]
        u_name = u_map[z_val]
        v_name = v_map[z_val]
        q_name = get(q_map, z_val, nothing)

        temp_series = _extract_time_series(ds, t_name, nt, station_index)
        u_series = _extract_time_series(ds, u_name, nt, station_index)
        v_series = _extract_time_series(ds, v_name, nt, station_index)
        q_series =
            q_name === nothing ? fill(DEFAULT_Q, nt) :
            _extract_time_series(ds, q_name, nt, station_index)

        finite_temps = filter(isfinite, temp_series)
        temp_k = temp_series .+ (!isempty(finite_temps) && maximum(finite_temps) < 200.0 ? 273.15 : 0.0)

        theta_series = temp_k .* (P0 ./ p_series) .^ R_OVER_CP

        @inbounds for it = 1:nt
            t_val = t_values[it]
            u_val = u_series[it]
            v_val = v_series[it]
            th_val = theta_series[it]
            q_val = q_series[it]
            if isfinite(t_val) && isfinite(u_val) && isfinite(v_val) && isfinite(th_val)
                push!(time_col, t_val)
                push!(z_col, z_val)
                push!(u_col, u_val)
                push!(v_col, v_val)
                push!(theta_col, th_val)
                push!(q_col, isfinite(q_val) ? q_val : DEFAULT_Q)
            end
        end
    end
end

function _as_time_height_matrix(values, nt::Int, nz::Int, varname::String)
    data = coalesce.(Array(values), NaN)
    if size(data) == (nt, nz)
        return Float64.(data)
    elseif size(data) == (nz, nt)
        return permutedims(Float64.(data), (2, 1))
    elseif ndims(data) == 1 && length(data) == nt
        return repeat(reshape(Float64.(data), nt, 1), 1, nz)
    elseif ndims(data) == 1 && length(data) == nz
        return repeat(reshape(Float64.(data), 1, nz), nt, 1)
    end

    throw(
        DimensionMismatch(
            "$varname has size $(size(data)); expected (time,height)=($nt,$nz) or (height,time)=($nz,$nt)",
        ),
    )
end

function _pressure_to_pa!(p::Matrix{Float64})
    if maximum(filter(isfinite, p); init=0.0) < 2_000.0
        p .*= 100.0
    end
    return p
end

function _build_theta_matrix(ds::NCDataset, nt::Int, nz::Int)
    theta_name = _resolve_var_name(ds, THETA_ALIASES; required=false)
    if theta_name !== nothing
        return _as_time_height_matrix(ds[theta_name][:], nt, nz, theta_name)
    end

    temp_name = _resolve_var_name(ds, TEMP_ALIASES)
    pressure_name = _resolve_var_name(ds, PRESSURE_ALIASES)

    t_raw = _as_time_height_matrix(ds[temp_name][:], nt, nz, temp_name)
    p_pa = _as_time_height_matrix(ds[pressure_name][:], nt, nz, pressure_name)
    _pressure_to_pa!(p_pa)

    finite_temps = filter(isfinite, t_raw)
    t_k = t_raw .+ (!isempty(finite_temps) && maximum(finite_temps) < 200.0 ? 273.15 : 0.0)

    return t_k .* (P0 ./ p_pa) .^ R_OVER_CP
end

function _build_wind_components(ds::NCDataset, nt::Int, nz::Int)
    u_name = _resolve_var_name(ds, U_ALIASES; required=false)
    v_name = _resolve_var_name(ds, V_ALIASES; required=false)

    if u_name !== nothing && v_name !== nothing
        u = _as_time_height_matrix(ds[u_name][:], nt, nz, u_name)
        v = _as_time_height_matrix(ds[v_name][:], nt, nz, v_name)
        return u, v
    end

    ws_name = _resolve_var_name(ds, WS_ALIASES)
    wd_name = _resolve_var_name(ds, WD_ALIASES)

    ws = _as_time_height_matrix(ds[ws_name][:], nt, nz, ws_name)
    wd = _as_time_height_matrix(ds[wd_name][:], nt, nz, wd_name)

    rad = deg2rad.(wd)
    u = -ws .* sin.(rad)
    v = -ws .* cos.(rad)
    return u, v
end

function _build_q_matrix(ds::NCDataset, nt::Int, nz::Int)
    q_name = _resolve_var_name(ds, Q_ALIASES; required=false)
    if q_name === nothing
        return fill(DEFAULT_Q, nt, nz)
    end
    q = _as_time_height_matrix(ds[q_name][:], nt, nz, q_name)
    q[.!isfinite.(q)] .= DEFAULT_Q
    return q
end

function convert_ncar_cases99_to_trajectory(
    input_dir::AbstractString,
    output_csv::AbstractString,
)
    isdir(input_dir) || throw(ArgumentError("Source directory not found: $input_dir"))

    nc_files = _find_netcdf_files(input_dir)
    isempty(nc_files) && throw(ArgumentError("No NetCDF files found in: $input_dir"))

    time_col = Float64[]
    z_col = Float64[]
    u_col = Float64[]
    v_col = Float64[]
    theta_col = Float64[]
    q_col = Float64[]

    station_index = try
        parse(Int, get(ENV, "CASES99_STATION_INDEX", "1"))
    catch
        1
    end

    for file in nc_files
        NCDataset(file, "r") do ds
            time_name = _resolve_var_name(ds, TIME_ALIASES)
            t_values = _read_time_seconds(ds, time_name)

            z_name = _resolve_var_name(ds, HEIGHT_ALIASES; required=false)
            if z_name !== nothing
                z_values = Float64.(collect(ds[z_name][:]))
                nt = length(t_values)
                nz = length(z_values)

                theta = _build_theta_matrix(ds, nt, nz)
                u, v = _build_wind_components(ds, nt, nz)
                q = _build_q_matrix(ds, nt, nz)

                @inbounds for it = 1:nt
                    for iz = 1:nz
                        t_val = t_values[it]
                        z_val = z_values[iz]
                        u_val = u[it, iz]
                        v_val = v[it, iz]
                        th_val = theta[it, iz]
                        q_val = q[it, iz]

                        if isfinite(t_val) &&
                           isfinite(z_val) &&
                           isfinite(u_val) &&
                           isfinite(v_val) &&
                           isfinite(th_val) &&
                           isfinite(q_val)
                            push!(time_col, t_val)
                            push!(z_col, z_val)
                            push!(u_col, u_val)
                            push!(v_col, v_val)
                            push!(theta_col, th_val)
                            push!(q_col, q_val)
                        end
                    end
                end
            else
                _build_rows_from_wide_schema!(
                    ds,
                    t_values,
                    station_index,
                    time_col,
                    z_col,
                    u_col,
                    v_col,
                    theta_col,
                    q_col,
                )
            end
        end
    end

    isempty(time_col) &&
        throw(ArgumentError("No valid rows were extracted from NetCDF sources"))

    t_min = minimum(time_col)
    time_col .-= t_min

    perm = sortperm(eachindex(time_col), by=i -> (time_col[i], z_col[i]))

    seen = Set{Tuple{Float64,Float64}}()
    unique_indices = Int[]
    sizehint!(unique_indices, length(perm))
    for idx in perm
        key = (time_col[idx], z_col[idx])
        if !(key in seen)
            push!(seen, key)
            push!(unique_indices, idx)
        end
    end

    z_counts = Dict{Float64,Int}()
    for idx in unique_indices
        t = time_col[idx]
        z_counts[t] = get(z_counts, t, 0) + 1
    end

    nz_required = maximum(values(z_counts))
    complete_times = Set{Float64}(t for (t, count) in z_counts if count == nz_required)

    valid_indices = filter(idx -> time_col[idx] in complete_times, unique_indices)

    mkpath(dirname(output_csv))
    CSV.write(
        output_csv,
        (
            time=time_col[valid_indices],
            z=z_col[valid_indices],
            u=u_col[valid_indices],
            v=v_col[valid_indices],
            theta=theta_col[valid_indices],
            q=q_col[valid_indices],
        ),
    )

    println("Generated CASES-99 trajectory CSV: $output_csv")
    println("  Rows: $(length(valid_indices))")
    println("  Unique times: $(length(complete_times))")
    println("  Unique z levels: $nz_required")
end

function main(argv::Vector{String})
    repo_root = normpath(joinpath(@__DIR__, ".."))

    input_dir = if !isempty(argv)
        argv[1]
    else
        get(ENV, "NCAR_CASES99_DIR", joinpath(repo_root, "data", "ncar_eol_dee0099881"))
    end

    output_csv = if length(argv) >= 2
        argv[2]
    else
        joinpath(repo_root, "data", "drafts", "trajectories", "trajectory_cases_99.csv")
    end

    if !isdir(input_dir)
        @error "Conversion failed: source directory not found: $input_dir"
        println("Set NCAR_CASES99_DIR or pass an explicit input directory.")
        exit(1)
    end

    convert_ncar_cases99_to_trajectory(input_dir, output_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end