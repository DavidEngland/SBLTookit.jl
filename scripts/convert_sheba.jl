#!/usr/bin/env julia
# scripts/convert_sheba.jl — Convert SHEBA Met City ASCII profiles to a standard trajectory CSV.
#
# Usage:
#   julia --project=. scripts/convert_sheba.jl [prof_file_all6_ed_hd.txt] [output_csv]
#
# Schema reference: prof_file_all6_ed_hd.txt provides 5-level tower profiles with
# JD, z1..z5, ws1..ws5, wd1..wd5, T1..T5, q1..q5, u*1..u*5, hs1..hs5, fl1..fl5.

using CSV
using Dates

const DEFAULT_Q = 0.005
const P0 = 100000.0
const R_OVER_CP = 0.286
const NUM_LEVELS = 5

# SHEBA Julian Day epoch: 1997-01-01T00:00:00 UTC = Day 1.0
const SHEBA_EPOCH_UNIX = datetime2unix(DateTime(1997, 1, 1, 0, 0, 0))

sheba_jd_to_unix(jd::Real) = SHEBA_EPOCH_UNIX + (jd - 1.0) * 86_400.0

function _level_column(prefix::AbstractString, level::Int)
    name = "$prefix$level"
    return Symbol(name)
end

function convert_sheba_profile_to_trajectory(
    input_txt::AbstractString,
    output_csv::AbstractString,
)
    isfile(input_txt) || throw(ArgumentError("Source file not found: $input_txt"))

    table = CSV.File(input_txt; delim = ' ', ignorerepeated = true, header = 1)

    time_col = Float64[]
    z_col = Float64[]
    u_col = Float64[]
    v_col = Float64[]
    theta_col = Float64[]
    q_col = Float64[]

    for row in table
        jd = getproperty(row, :JD)
        isfinite(jd) || continue
        t_val = sheba_jd_to_unix(jd)

        for level = 1:NUM_LEVELS
            fl = getproperty(row, _level_column("fl", level))
            # QC failure: sonic/turbulence flag set means this level is unreliable.
            fl == 1 && continue

            z_val = Float64(getproperty(row, _level_column("z", level)))
            ws_val = Float64(getproperty(row, _level_column("ws", level)))
            wd_val = Float64(getproperty(row, _level_column("wd", level)))
            temp_c = Float64(getproperty(row, _level_column("T", level)))
            q_gkg = Float64(getproperty(row, _level_column("q", level)))

            temp_k = temp_c + 273.15
            q_val = q_gkg * 1.0e-3
            theta_val = temp_k * (P0 / P0)^R_OVER_CP

            rad = deg2rad(wd_val)
            u_val = -ws_val * sin(rad)
            v_val = -ws_val * cos(rad)

            if isfinite(t_val) &&
               isfinite(z_val) &&
               isfinite(u_val) &&
               isfinite(v_val) &&
               isfinite(theta_val)
                push!(time_col, t_val)
                push!(z_col, z_val)
                push!(u_col, u_val)
                push!(v_col, v_val)
                push!(theta_col, theta_val)
                push!(q_col, isfinite(q_val) ? q_val : DEFAULT_Q)
            end
        end
    end

    isempty(time_col) &&
        throw(ArgumentError("No valid rows were extracted from SHEBA profile source"))

    t_min = minimum(time_col)
    time_col .-= t_min

    perm = sortperm(eachindex(time_col), by = i -> (time_col[i], z_col[i]))

    # Deduplicate (time, z) pairs in sorted order without intermediate vectors.
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

    # Count vertical levels present per timestamp to find complete profiles.
    z_counts = Dict{Float64,Int}()
    for idx in unique_indices
        t = time_col[idx]
        z_counts[t] = get(z_counts, t, 0) + 1
    end

    nz_required = length(unique(z_col[unique_indices]))
    complete_times = Set{Float64}(t for (t, count) in z_counts if count == nz_required)

    # Keep only complete timestamps so the benchmark loader receives a full
    # rectangular (time, z) grid with no missing state entries.
    valid_indices = filter(idx -> time_col[idx] in complete_times, unique_indices)

    mkpath(dirname(output_csv))
    CSV.write(
        output_csv,
        (
            time = time_col[valid_indices],
            z = z_col[valid_indices],
            u = u_col[valid_indices],
            v = v_col[valid_indices],
            theta = theta_col[valid_indices],
            q = q_col[valid_indices],
        ),
    )

    println("Generated SHEBA trajectory CSV: $output_csv")
    println("  Rows: $(length(valid_indices))")
    println("  Unique times: $(length(complete_times))")
    println("  Unique z levels: $nz_required")
end

function main(argv::Vector{String})
    repo_root = normpath(joinpath(@__DIR__, ".."))

    input_txt = if !isempty(argv)
        argv[1]
    else
        get(
            ENV,
            "SHEBA_PROFILE_FILE",
            joinpath(repo_root, "data", "raw", "sheba", "prof_file_all6_ed_hd.txt"),
        )
    end

    output_csv = if length(argv) >= 2
        argv[2]
    else
        joinpath(repo_root, "data", "drafts", "trajectories", "trajectory_sheba.csv")
    end

    if !isfile(input_txt)
        println("Conversion skipped: source file not found: $input_txt")
        println("Set SHEBA_PROFILE_FILE or pass an explicit input file path.")
        return
    end

    convert_sheba_profile_to_trajectory(input_txt, output_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
