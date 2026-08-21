#!/usr/bin/env julia
# scripts/convert_sheba.jl — Convert SHEBA Met City ASCII profiles to a standard trajectory CSV.

using CSV
using Dates
using Logging

const DEFAULT_Q = 0.005
const P0 = 100_000.0        # Reference pressure (Pa)
# const R_OVER_CP = 0.2854     # R_d / c_p for dry air
const G_ACCEL = 9.80665      # Gravitational acceleration (m/s^2)
const CP_AIR = 1004.67       # Specific heat capacity at constant pressure (J/kg/K)
const NUM_LEVELS = 5

# SHEBA Julian Day epoch: 1997-01-01T00:00:00 UTC = Day 1.0
const SHEBA_EPOCH_UNIX = datetime2unix(DateTime(1997, 1, 1, 0, 0, 0))

sheba_jd_to_unix(jd::Real) = SHEBA_EPOCH_UNIX + (jd - 1.0) * 86_400.0

tofloat(x) = x === missing ? NaN : Float64(x)

function _level_column(prefix::AbstractString, level::Int)
    return Symbol("$prefix$level")
end

function convert_sheba_profile_to_trajectory(
    input_txt::AbstractString,
    output_csv::AbstractString,
)
    isfile(input_txt) || throw(ArgumentError("Source file not found: $input_txt"))

    table = CSV.File(
        input_txt;
        delim=('\t'),
        ignorerepeated=true,
        header=1,
        skipto=3,
        missingstring=["999", "9999", "999.0", "9999.0"],
    )

    time_col = Float64[]
    z_col = Float64[]
    u_col = Float64[]
    v_col = Float64[]
    theta_col = Float64[]
    q_col = Float64[]

    for row in table
        jd = tofloat(getproperty(row, :JD))
        isfinite(jd) || continue
        t_val = sheba_jd_to_unix(jd)

        for level = 1:NUM_LEVELS
            fl = tofloat(getproperty(row, _level_column("fl", level)))
            # Rejection: sonic/turbulence flag set OR flag value missing (NaN)
            (isnan(fl) || fl == 1) && continue

            z_val = tofloat(getproperty(row, _level_column("z", level)))
            ws_val = tofloat(getproperty(row, _level_column("ws", level)))
            wd_val = tofloat(getproperty(row, _level_column("wd", level)))
            temp_c = tofloat(getproperty(row, _level_column("T", level)))
            q_gkg = tofloat(getproperty(row, _level_column("q", level)))

            temp_k = temp_c + 273.15
            q_val = q_gkg * 1.0e-3

            # Hydrostatic pressure correction: theta = T * exp(g * z / (c_p * T))
            theta_val = temp_k * exp((G_ACCEL * z_val) / (CP_AIR * temp_k))

            rad = deg2rad(wd_val)
            u_val = -ws_val * sin(rad)
            v_val = -ws_val * cos(rad)

            if isfinite(z_val) &&
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

    nz_required = NUM_LEVELS
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
        @error "Conversion failed: source file not found: $input_txt"
        println("Set SHEBA_PROFILE_FILE or pass an explicit input file path.")
        exit(1)
    end

    convert_sheba_profile_to_trajectory(input_txt, output_csv)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end