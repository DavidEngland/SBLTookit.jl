# src/GSPTTimeLoop.jl
module GSPTTimeLoop

using NCDatasets, DynamicQuantities, Plots, Dates, Statistics
using ..GSPTPhase2

export ingest_netcdf_gspt, plot_gspt_transition

"""
    get_nc_var(ds::NCDataset, candidate_names::Vector{String}; required=true)

Extracts variables from NetCDF datasets safely without forcing DateTimes to Float64.
"""
function get_nc_var(ds::NCDataset, candidate_names::Vector{String}; required=true)
    for key in candidate_names
        if haskey(ds, key)
            val = ds[key][:]
            if val isa AbstractArray{<:Number}
                return Float64.(val)
            elseif val isa AbstractArray{<:Union{Missing,Number}}
                return Float64.(coalesce.(val, NaN))
            else
                return val
            end
        end
    end
    if required
        available = join(keys(ds), ", ")
        error("None of the keys $(candidate_names) found in NetCDF dataset. Available keys: [$available]")
    end
    return nothing
end

"""
    ensure_2d_matrix(arr, N_z::Int, N_t::Int)

Normalizes 1D vectors and transposed 2D arrays into a standard (N_z, N_t) matrix.
"""
function ensure_2d_matrix(arr, N_z::Int, N_t::Int)
    arr === nothing && return zeros(N_z, N_t)

    if ndims(arr) == 1
        vec_arr = vec(arr)
        mat = zeros(N_z, N_t)
        if length(vec_arr) == N_t
            mat[1, :] .= vec_arr  # Place surface time-series in row 1
        elseif length(vec_arr) == N_z
            mat .= repeat(reshape(vec_arr, N_z, 1), 1, N_t)
        end
        return mat
    elseif ndims(arr) == 2
        s1, s2 = size(arr)
        if s1 == N_z && s2 == N_t
            return arr
        elseif s1 == N_t && s2 == N_z
            return permutedims(arr, (2, 1))
        else
            mat = zeros(N_z, N_t)
            min_z = min(s1, N_z)
            min_t = min(s2, N_t)
            mat[1:min_z, 1:min_t] .= arr[1:min_z, 1:min_t]
            return mat
        end
    end
    return zeros(N_z, N_t)
end

"""
    try_extract_tower_2d(ds::NCDataset)

Parses level-suffixed tower observations (e.g. u_5m, u_10m, tc_1_5m) into 2D (z, t) arrays.
"""
function try_extract_tower_2d(ds::NCDataset)
    var_names = keys(ds)
    suffix_map = Dict{Float64,String}()

    for name in var_names
        m = match(r"_(?:u|v|w|tc|T|U|V|w_tc|u_w|v_w)_?(\d+(?:[\._]\d+)?)m$", name)
        if m !== nothing
            h_str = m.captures[1]
            val = parse(Float64, replace(h_str, "_" => "."))
            suffix_map[val] = "_" * h_str * "m"
        end
    end

    isempty(suffix_map) && return nothing

    sorted_z = sort(collect(keys(suffix_map)))
    N_z = length(sorted_z)

    time_raw = get_nc_var(ds, ["time", "t", "base_time"])
    N_t = length(time_raw)

    u_mat = zeros(N_z, N_t)
    v_mat = zeros(N_z, N_t)
    th_mat = zeros(N_z, N_t)
    wth_mat = zeros(N_z, N_t)
    uw_mat = zeros(N_z, N_t)
    vw_mat = zeros(N_z, N_t)

    for (i, z_val) in enumerate(sorted_z)
        sfx = suffix_map[z_val]

        get_level_data = (candidates) -> begin
            for c in candidates
                for key_cand in [c * sfx, c * replace(sfx, r"^_" => "")]
                    if haskey(ds, key_cand)
                        d = ds[key_cand][:]
                        vec_d = d isa AbstractArray{<:Union{Missing,Number}} ? Float64.(vec(coalesce.(d, NaN))) : Float64.(vec(d))

                        if length(vec_d) == N_t
                            return vec_d
                        elseif length(vec_d) > N_t
                            return vec_d[1:N_t]
                        elseif length(vec_d) > 0
                            res = zeros(N_t)
                            res[1:length(vec_d)] .= vec_d
                            return res
                        end
                    end
                end
            end
            return zeros(N_t)
        end

        u_arr = get_level_data(["u", "U", "u_wind"])
        v_arr = get_level_data(["v", "V", "v_wind"])
        tc_arr = get_level_data(["tc", "T", "temp", "t"])
        wth_arr = get_level_data(["w_tc", "wthv", "wth", "wt"])
        uw_arr = get_level_data(["u_w", "uw"])
        vw_arr = get_level_data(["v_w", "vw"])

        valid_tc = filter(!isnan, tc_arr)
        th_arr = (!isempty(valid_tc) && mean(valid_tc) < 100.0) ? tc_arr .+ 273.15 : tc_arr

        u_mat[i, :] .= u_arr
        v_mat[i, :] .= v_arr
        th_mat[i, :] .= th_arr
        wth_mat[i, :] .= wth_arr
        uw_mat[i, :] .= uw_arr
        vw_mat[i, :] .= vw_arr
    end

    return (z=sorted_z, time=time_raw, u=u_mat, v=v_mat, th=th_mat, wth=wth_mat, uw=uw_mat, vw=vw_mat)
end

"""
    ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
"""
function ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
    NCDataset(nc_path, "r") do ds
        z = get_nc_var(ds, ["zf", "zt", "z", "height", "heights", "lev", "level", "alt", "z_m"]; required=false)

        local time_raw, u_mat, v_mat, th_mat, wth_mat, uw_mat, vw_mat, tke_mat, Km_mat

        if z === nothing
            tower_data = try_extract_tower_2d(ds)
            if tower_data === nothing
                available = join(keys(ds), ", ")
                error("No vertical coordinate found and could not parse level suffixes. Available keys: [$available]")
            end
            z = tower_data.z
            time_raw = tower_data.time
            u_mat = tower_data.u
            v_mat = tower_data.v
            th_mat = tower_data.th
            wth_mat = tower_data.wth
            uw_mat = tower_data.uw
            vw_mat = tower_data.vw
            tke_mat = zeros(length(z), length(time_raw))
            Km_mat = zeros(length(z), length(time_raw))
        else
            time_raw = get_nc_var(ds, ["time", "t", "datetime"])
            N_z, N_t = length(z), length(time_raw)

            u_raw = get_nc_var(ds, ["u", "U", "u_wind", "wind_u", "eastward_wind"])
            v_raw = get_nc_var(ds, ["v", "V", "v_wind", "wind_v", "northward_wind"])
            th_raw = get_nc_var(ds, ["theta_v", "theta", "th", "THETA", "thv", "temp", "t"])
            wth_raw = get_nc_var(ds, ["wthv", "wth", "w_theta", "w_th", "flux_th", "wt", "heat_flux", "shf"])
            uw_raw = get_nc_var(ds, ["uw", "u_w", "uw_flux", "momentum_flux_u"])
            vw_raw = get_nc_var(ds, ["vw", "v_w", "vw_flux", "momentum_flux_v"])

            tke_raw = get_nc_var(ds, ["tke", "TKE", "q2"]; required=false)
            Km_raw = get_nc_var(ds, ["Km", "km", "KM", "K_m"]; required=false)

            u_mat = ensure_2d_matrix(u_raw, N_z, N_t)
            v_mat = ensure_2d_matrix(v_raw, N_z, N_t)
            th_mat = ensure_2d_matrix(th_raw, N_z, N_t)
            wth_mat = ensure_2d_matrix(wth_raw, N_z, N_t)
            uw_mat = ensure_2d_matrix(uw_raw, N_z, N_t)
            vw_mat = ensure_2d_matrix(vw_raw, N_z, N_t)
            tke_mat = ensure_2d_matrix(tke_raw, N_z, N_t)
            Km_mat = ensure_2d_matrix(Km_raw, N_z, N_t)

            # Convert surface sensible heat flux (W/m^2) to kinematic flux (K m/s) if applicable
            if wth_raw !== nothing && maximum(abs.(wth_mat)) > 5.0
                wth_mat .= wth_mat ./ 1200.0
            end

            # Convert Celsius to Kelvin if needed
            valid_th = filter(!isnan, th_mat)
            if !isempty(valid_th) && mean(valid_th) < 100.0
                th_mat .+= 273.15
            end
        end

        N_z = length(z)
        N_t = length(time_raw)

        R_coord_2d = fill(NaN, N_z, N_t)
        C_const_2d = fill(NaN, N_z, N_t)
        C_coord_2d = fill(NaN, N_z, N_t)
        delta_clos_2d = fill(NaN, N_z, N_t)
        mask_2d = fill(false, N_z, N_t)

        for t in 1:N_t
            # Nocturnal filter: negative surface kinematic heat flux
            is_nocturnal = wth_mat[1, t] < 0.0
            !is_nocturnal && continue

            data_t = ProfileData(
                z,
                u_mat[:, t],
                v_mat[:, t],
                th_mat[:, t],
                wth_mat[:, t],
                uw_mat[:, t],
                vw_mat[:, t],
                tke_mat[:, t],
                Km_mat[:, t],
                σ_u, σ_v, σ_th, 0.01
            )

            res_t = compute_gspt(data_t; is_observation=true, S2_min=S2_min)

            R_coord_2d[:, t] .= res_t.const_geom.R_coord
            C_const_2d[:, t] .= res_t.const_geom.C_const
            C_coord_2d[:, t] .= res_t.const_geom.C_coord
            delta_clos_2d[:, t] .= res_t.obs_diag.closure_residual
            mask_2d[:, t] .= res_t.obs_diag.ill_conditioned_mask
        end

        return (z=z, time=time_raw, R_coord=R_coord_2d, C_const=C_const_2d,
            C_coord=C_coord_2d, delta_closure=delta_clos_2d, mask=mask_2d)
    end
end

"""
    plot_gspt_transition(results; save_path="gspt_transition_heatmap.png")
"""
function plot_gspt_transition(results; save_path="gspt_transition_heatmap.png")
    z = results.z
    t = results.time
    R_mat = results.R_coord

    p = heatmap(
        1:length(t), z, R_mat,
        clims=(-2.0, 2.0),
        color=:coolwarm,
        xlabel="Time Index / Step",
        ylabel="Height z (m)",
        title="Dynamic Nocturnal GSPT Transition Surface (R_coord)",
        colorbar_title="R_coord",
        size=(1000, 500)
    )

    savefig(p, save_path)
    return p
end

end # module GSPTTimeLoop