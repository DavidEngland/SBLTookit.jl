module GSPTTimeLoop

using NCDatasets, Plots, Dates, Statistics
using ..GSPTPhase2

export ingest_netcdf_gspt, plot_gspt_transition

# Atmospheric thermodynamic constants
const RHO_AIR = 1.225     # Reference air density [kg/m^3]
const CP_AIR = 1004.67   # Specific heat of dry air [J/(kg K)]
const RHO_CP = RHO_AIR * CP_AIR  # ~1230.7 J/(m^3 K)
const LAPSE_DRY = 0.0098  # Dry adiabatic lapse rate [K/m]

# Match level-suffixed tower variables explicitly: variable + height + optional 'm'
const TOWER_LEVEL_RE = r"^(w_tc|wthv|wth|wt|u_w|v_w|hsb|usb|hlb|ws|wd|tc|temp|theta_v|theta|th|t|u|v|w|T|U|V)_?(\d+(?:[\._]\d+)?)m?$"

# ---------------------------------------------------------------------------
# Layer 1: Data Extraction & Variable Identification
# ---------------------------------------------------------------------------

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
        error("None of candidate keys $(candidate_names) found in dataset. Available keys: [$available]")
    end
    return nothing
end

"""
    extract_decoded_time(ds::NCDataset)

Extracts and decodes NetCDF time variables into `Vector{DateTime}` or `Vector{Float64}`.
"""
function extract_decoded_time(ds::NCDataset)
    time_raw = get_nc_var(ds, ["time", "t", "datetime", "base_time"])
    if time_raw isa Vector{<:Dates.AbstractDateTime}
        return time_raw
    end

    # Try attribute decoding if raw numeric values returned
    if haskey(ds, "time") && haskey(ds["time"].attrib, "units")
        try
            return NCDatasets.num2date(time_raw, ds["time"].attrib["units"])
        catch
            return time_raw
        end
    end
    return time_raw
end

"""
    ensure_2d_matrix(arr, N_z::Int, N_t::Int; name::String="variable")

Initializes missing array elements with NaN to avoid manufacturing false zero gradients.
"""
function ensure_2d_matrix(arr, N_z::Int, N_t::Int; name::String="variable")
    mat = fill(NaN, N_z, N_t)
    arr === nothing && return mat

    if ndims(arr) == 1
        vec_arr = vec(arr)
        if length(vec_arr) == N_t
            mat[1, :] .= vec_arr  # Surface time series in row 1
        elseif length(vec_arr) == N_z
            mat .= repeat(reshape(vec_arr, N_z, 1), 1, N_t)
        else
            @warn "ensure_2d_matrix: 1D '$name' has length $(length(vec_arr)), matching neither N_z=$N_z nor N_t=$N_t; filling with NaN."
        end
    elseif ndims(arr) == 2
        s1, s2 = size(arr)
        if s1 == N_z && s2 == N_t
            return arr
        elseif s1 == N_t && s2 == N_z
            return permutedims(arr, (2, 1))
        else
            @warn "ensure_2d_matrix: 2D '$name' has shape $(size(arr)), matching neither ($N_z,$N_t) nor ($N_t,$N_z); truncating/padding with NaN."
            min_z = min(s1, N_z)
            min_t = min(s2, N_t)
            mat[1:min_z, 1:min_t] .= arr[1:min_z, 1:min_t]
        end
    end
    return mat
end

# ---------------------------------------------------------------------------
# Layer 2: Physical & Thermodynamic Normalization
# ---------------------------------------------------------------------------

"""
    normalize_thermodynamics(ds::NCDataset, z::Vector{Float64}, N_t::Int)

Extracts and converts thermodynamic variables to Virtual Potential Temperature (\\theta_v).
Calculates \\theta_v = \\theta (1 + 0.61 q) or \\theta_v \\approx T + \\gamma_{dry} z as fallback.
"""
function normalize_thermodynamics(ds::NCDataset, z::Vector{Float64}, N_t::Int)
    N_z = length(z)

    # 1. Direct Virtual Potential Temperature (\theta_v)
    thv_raw = get_nc_var(ds, ["theta_v", "thv", "THETA_V"]; required=false)
    if thv_raw !== nothing
        thv_mat = ensure_2d_matrix(thv_raw, N_z, N_t; name="theta_v")
        valid = filter(!isnan, thv_mat)
        if !isempty(valid) && mean(valid) < 100.0
            thv_mat .+= 273.15
        end
        return thv_mat
    end

    # 2. Potential Temperature (\theta)
    th_raw = get_nc_var(ds, ["theta", "th", "THETA"]; required=false)
    q_raw = get_nc_var(ds, ["q", "qv", "SPEC_HUM"]; required=false)

    if th_raw !== nothing
        th_mat = ensure_2d_matrix(th_raw, N_z, N_t; name="theta")
        valid = filter(!isnan, th_mat)
        if !isempty(valid) && mean(valid) < 100.0
            th_mat .+= 273.15
        end

        if q_raw !== nothing
            q_mat = ensure_2d_matrix(q_raw, N_z, N_t; name="q")
            return th_mat .* (1.0 .+ 0.61 .* q_mat)
        end
        return th_mat
    end

    # 3. Absolute Temperature (T) conversion to \theta_v
    t_raw = get_nc_var(ds, ["temp", "t", "T", "tc"]; required=true)
    t_mat = ensure_2d_matrix(t_raw, N_z, N_t; name="temperature")
    valid = filter(!isnan, t_mat)
    if !isempty(valid) && mean(valid) < 100.0
        t_mat .+= 273.15
    end

    thv_mat = fill(NaN, N_z, N_t)
    for i in 1:N_z
        thv_mat[i, :] .= t_mat[i, :] .+ (LAPSE_DRY * z[i])
    end
    return thv_mat
end

"""
    normalize_heat_flux(wth_mat; flux_convention=:auto, rho_cp=RHO_CP)

Normalizes heat flux to kinematic units (K m/s). Converts Sensible Heat Flux (W/m^2) using \\rho c_p.
"""
function normalize_heat_flux(wth_mat::Matrix{Float64}; flux_convention::Symbol=:auto, rho_cp::Float64=RHO_CP)
    wth_norm = copy(wth_mat)
    valid_vals = filter(!isnan, wth_norm)
    isempty(valid_vals) && return wth_norm

    is_sensible = flux_convention == :sensible || (flux_convention == :auto && maximum(abs.(valid_vals)) > 5.0)
    if is_sensible
        wth_norm ./= rho_cp
    end
    return wth_norm
end

"""
    normalize_momentum_flux(uw_raw, vw_raw, N_z::Int, N_t::Int)

Converts friction velocity u_* (usb) elementwise via \\overline{u'w'} = -u_*^2 for valid positive values.
"""
function normalize_momentum_flux(ds::NCDataset, N_z::Int, N_t::Int)
    uw_raw = get_nc_var(ds, ["uw", "u_w", "uw_flux", "momentum_flux_u"]; required=false)
    vw_raw = get_nc_var(ds, ["vw", "v_w", "vw_flux", "momentum_flux_v"]; required=false)
    usb_raw = get_nc_var(ds, ["usb", "ustar", "u_star"]; required=false)

    uw_mat = ensure_2d_matrix(uw_raw, N_z, N_t; name="uw")
    vw_mat = ensure_2d_matrix(vw_raw, N_z, N_t; name="vw")

    if usb_raw !== nothing
        usb_mat = ensure_2d_matrix(usb_raw, N_z, N_t; name="usb")
        for t in 1:N_t, z_idx in 1:N_z
            val = usb_mat[z_idx, t]
            if isnan(uw_mat[z_idx, t]) && !isnan(val) && val >= 0.0
                uw_mat[z_idx, t] = -(val^2)
            end
        end
    end
    return uw_mat, vw_mat
end

# ---------------------------------------------------------------------------
# Tower & Gridded Ingestion Workflows
# ---------------------------------------------------------------------------

function try_extract_tower_2d(ds::NCDataset)
    var_names = keys(ds)
    suffix_map = Dict{Float64,String}()

    for name in var_names
        m = match(TOWER_LEVEL_RE, name)
        if m !== nothing
            h_str = m.captures[2]
            val = parse(Float64, replace(h_str, "_" => "."))
            suffix_map[val] = h_str
        end
    end

    isempty(suffix_map) && return nothing

    sorted_z = sort(collect(keys(suffix_map)))
    N_z = length(sorted_z)
    time_raw = extract_decoded_time(ds)
    N_t = length(time_raw)

    u_mat = fill(NaN, N_z, N_t)
    v_mat = fill(NaN, N_z, N_t)
    tc_mat = fill(NaN, N_z, N_t)
    wth_mat = fill(NaN, N_z, N_t)
    uw_mat = fill(NaN, N_z, N_t)
    vw_mat = fill(NaN, N_z, N_t)

    for (i, z_val) in enumerate(sorted_z)
        h_str = suffix_map[z_val]
        sfx_candidates = ["_" * h_str * "m", h_str * "m", "_" * h_str, h_str]

        get_level_data = (candidates) -> begin
            for c in candidates
                for sfx in sfx_candidates
                    key_cand = c * sfx
                    if haskey(ds, key_cand)
                        d = ds[key_cand][:]
                        vec_d = d isa AbstractArray{<:Union{Missing,Number}} ? Float64.(vec(coalesce.(d, NaN))) : Float64.(vec(d))
                        return length(vec_d) >= N_t ? vec_d[1:N_t] : [vec_d; fill(NaN, N_t - length(vec_d))]
                    end
                end
            end
            return fill(NaN, N_t)
        end

        u_mat[i, :] .= get_level_data(["u", "U"])
        v_mat[i, :] .= get_level_data(["v", "V"])
        tc_mat[i, :] .= get_level_data(["tc", "temp", "t", "T"])
        wth_mat[i, :] .= get_level_data(["w_tc", "wthv", "wth", "wt", "hsb"])
        uw_mat[i, :] .= get_level_data(["u_w", "uw"])
        vw_mat[i, :] .= get_level_data(["v_w", "vw"])

        # Elementwise friction velocity fallback
        usb_arr = get_level_data(["usb", "ustar"])
        for t in 1:N_t
            if isnan(uw_mat[i, t]) && !isnan(usb_arr[t]) && usb_arr[t] >= 0.0
                uw_mat[i, t] = -(usb_arr[t]^2)
            end
        end
    end

    # Thermodynamic conversion T -> \theta_v
    valid_t = filter(!isnan, tc_mat)
    if !isempty(valid_t) && mean(valid_t) < 100.0
        tc_mat .+= 273.15
    end
    thv_mat = fill(NaN, N_z, N_t)
    for i in 1:N_z
        thv_mat[i, :] .= tc_mat[i, :] .+ (LAPSE_DRY * sorted_z[i])
    end

    return (z=sorted_z, time=time_raw, u=u_mat, v=v_mat, thv=thv_mat, wth=wth_mat, uw=uw_mat, vw=vw_mat)
end

# ---------------------------------------------------------------------------
# Main Entry Point & Validation Engine
# ---------------------------------------------------------------------------

"""
    ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05,
                        flux_convention=:auto, rho_cp=RHO_CP, nocturnal_only=true,
                        use_threads=false)

Ingests NetCDF observations, normalizes thermodynamics and physical fluxes, performs profile validation,
and executes GSPT over valid time columns. Returns explicit diagnostic masks.
"""
function ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05,
    flux_convention::Symbol=:auto, rho_cp::Float64=RHO_CP, nocturnal_only::Bool=true,
    use_threads::Bool=false)

    NCDataset(nc_path, "r") do ds
        z_raw = get_nc_var(ds, ["zf", "zt", "z", "height", "heights", "lev", "level", "alt", "z_m"]; required=false)

        local time_raw, u_mat, v_mat, thv_mat, wth_mat, uw_mat, vw_mat, tke_mat, Km_mat, z

        if z_raw === nothing
            tower_data = try_extract_tower_2d(ds)
            tower_data === nothing && error("Could not parse vertical coordinate or tower level suffixes in $(nc_path).")
            z = tower_data.z
            time_raw = tower_data.time
            u_mat = tower_data.u
            v_mat = tower_data.v
            thv_mat = tower_data.thv
            wth_mat = tower_data.wth
            uw_mat = tower_data.uw
            vw_mat = tower_data.vw
            tke_mat = fill(NaN, length(z), length(time_raw))
            Km_mat = fill(NaN, length(z), length(time_raw))
        else
            time_raw = extract_decoded_time(ds)
            N_t = length(time_raw)

            # Extract 1D height coordinate along vertical dimension
            if ndims(z_raw) == 2
                z = size(z_raw, 2) == N_t ? Float64.(z_raw[:, 1]) : Float64.(z_raw[1, :])
            else
                z = Float64.(vec(z_raw))
            end
            N_z = length(z)

            u_raw = get_nc_var(ds, ["u", "U", "u_wind", "eastward_wind"])
            v_raw = get_nc_var(ds, ["v", "V", "v_wind", "northward_wind"])
            wth_raw = get_nc_var(ds, ["wthv", "wth", "w_theta", "w_th", "flux_th", "wt", "heat_flux", "shf", "hsb"])
            tke_raw = get_nc_var(ds, ["tke", "TKE", "q2"]; required=false)
            Km_raw = get_nc_var(ds, ["Km", "km", "KM", "K_m"]; required=false)

            u_mat = ensure_2d_matrix(u_raw, N_z, N_t; name="u")
            v_mat = ensure_2d_matrix(v_raw, N_z, N_t; name="v")
            wth_mat = ensure_2d_matrix(wth_raw, N_z, N_t; name="wth")
            tke_mat = ensure_2d_matrix(tke_raw, N_z, N_t; name="tke")
            Km_mat = ensure_2d_matrix(Km_raw, N_z, N_t; name="Km")

            thv_mat = normalize_thermodynamics(ds, z, N_t)
            uw_mat, vw_mat = normalize_momentum_flux(ds, N_z, N_t)
        end

        # Normalize heat flux to kinematic units (K m/s)
        wth_mat = normalize_heat_flux(wth_mat; flux_convention=flux_convention, rho_cp=rho_cp)

        # Ensure monotonic height ordering
        if !issorted(z)
            p = sortperm(z)
            z = z[p]
            u_mat = u_mat[p, :];
            v_mat = v_mat[p, :];
            thv_mat = thv_mat[p, :]
            wth_mat = wth_mat[p, :];
            uw_mat = uw_mat[p, :];
            vw_mat = vw_mat[p, :]
            tke_mat = tke_mat[p, :];
            Km_mat = Km_mat[p, :]
        end

        N_z, N_t = length(z), length(time_raw)

        R_coord_2d = fill(NaN, N_z, N_t)
        C_const_2d = fill(NaN, N_z, N_t)
        C_coord_2d = fill(NaN, N_z, N_t)
        delta_clos_2d = fill(NaN, N_z, N_t)

        # Explicit Diagnostic Masks
        mask_missing = fill(false, N_z, N_t)
        mask_nocturnal = fill(false, N_z, N_t)
        mask_ill_conditioned = fill(false, N_z, N_t)

        compute_column! = t -> begin
            sfc_flux = wth_mat[1, t]

            # Nocturnal Masking (\overline{w'\theta'} < 0)
            is_noct = !isnan(sfc_flux) && sfc_flux < 0.0
            mask_nocturnal[:, t] .= is_noct
            if nocturnal_only && !is_noct
                return nothing
            end

            # Check profile completeness (reject NaNs in critical fields)
            u_col, v_col, thv_col = u_mat[:, t], v_mat[:, t], thv_mat[:, t]
            has_missing = any(isnan, u_col) || any(isnan, v_col) || any(isnan, thv_col)
            if has_missing
                mask_missing[:, t] .= true
                return nothing
            end

            data_t = ProfileData(
                z, u_col, v_col, thv_col,
                wth_mat[:, t], uw_mat[:, t], vw_mat[:, t],
                tke_mat[:, t], Km_mat[:, t],
                σ_u, σ_v, σ_th, 0.01
            )
            # 1. Fix 2D/N-D Vertical Coordinate Dimensionality
            if ndims(z_raw) >= 2
                s1, s2 = size(z_raw)[1:2]
                if s2 == N_t
                    z = Float64.(z_raw[:, 1])
                elseif s1 == N_t
                    z = Float64.(z_raw[1, :])
                else
                    z = Float64.(z_raw[:, 1])
                end
            else
                z = Float64.(z_raw)
            end
            N_z = length(z)

            # 2. Add Reshaping Support for Multi-Level 1D Flux Arrays (e.g., 5 levels * 144 times = 720)
            function ensure_2d_matrix(arr, N_z::Int, N_t::Int; name::String="variable")
                mat = fill(NaN, N_z, N_t)
                arr === nothing && return mat

                if ndims(arr) == 1
                    vec_arr = vec(arr)
                    if length(vec_arr) == N_t
                        mat[1, :] .= vec_arr
                    elseif length(vec_arr) == N_z
                        mat .= repeat(reshape(vec_arr, N_z, 1), 1, N_t)
                    elseif length(vec_arr) == N_z * N_t
                        return reshape(Float64.(vec_arr), N_z, N_t)
                    elseif mod(length(vec_arr), N_t) == 0
                        n_sub = div(length(vec_arr), N_t)
                        sub_mat = reshape(Float64.(vec_arr), n_sub, N_t)
                        min_z = min(n_sub, N_z)
                        mat[1:min_z, :] .= sub_mat[1:min_z, :]
                    else
                        @warn "ensure_2d_matrix: 1D '$name' length $(length(vec_arr)) unresolvable; filling with NaN."
                    end
                elseif ndims(arr) == 2
                    s1, s2 = size(arr)
                    if s1 == N_z && s2 == N_t
                        return arr
                    elseif s1 == N_t && s2 == N_z
                        return permutedims(arr, (2, 1))
                    end
                end
                return mat
            end

            # 3. Guard GSPT Derivative Engine Against Under-Determined Profiles (N_z < 3)
            if N_z < 3
                @warn "Dataset '$nc_path' has N_z = $N_z (< 3 levels). GSPT finite differences require N_z >= 3. Skipping solver execution."
                return (
                    z=z, time=time_raw,
                    R_coord=fill(NaN, N_z, N_t), C_const=fill(NaN, N_z, N_t),
                    C_coord=fill(NaN, N_z, N_t), delta_closure=fill(NaN, N_z, N_t),
                    mask_missing=fill(true, N_z, N_t), mask_nocturnal=fill(false, N_z, N_t),
                    mask_ill_conditioned=fill(false, N_z, N_t)
                )
            end
            res_t = compute_gspt(data_t; is_observation=true, S2_min=S2_min)

            R_coord_2d[:, t] .= res_t.const_geom.R_coord
            C_const_2d[:, t] .= res_t.const_geom.C_const
            C_coord_2d[:, t] .= res_t.const_geom.C_coord
            delta_clos_2d[:, t] .= res_t.obs_diag.closure_residual
            mask_ill_conditioned[:, t] .= res_t.obs_diag.ill_conditioned_mask
            return nothing
        end

        if use_threads
            Threads.@threads for t in 1:N_t
                compute_column!(t)
            end
        else
            for t in 1:N_t
                compute_column!(t)
            end
        end

        return (
            z=z, time=time_raw, R_coord=R_coord_2d, C_const=C_const_2d,
            C_coord=C_coord_2d, delta_closure=delta_clos_2d,
            mask_missing=mask_missing, mask_nocturnal=mask_nocturnal,
            mask_ill_conditioned=mask_ill_conditioned
        )
    end
end

# ---------------------------------------------------------------------------
# Layer 4: Transition Diagnostics & Visualization
# ---------------------------------------------------------------------------

"""
    plot_gspt_transition(results; save_path="gspt_transition_heatmap.png",
                          mask_ill_conditioned=true, mask_missing=true)

Plots dynamic R_coord heatmap over height and time. Masks ill-conditioned or missing observations.
"""
function plot_gspt_transition(results; save_path="gspt_transition_heatmap.png",
    mask_ill_conditioned::Bool=true, mask_missing::Bool=true)

    z = results.z
    t = results.time
    R_mat = copy(results.R_coord)

    if mask_ill_conditioned && hasproperty(results, :mask_ill_conditioned)
        R_mat[results.mask_ill_conditioned] .= NaN
    end
    if mask_missing && hasproperty(results, :mask_missing)
        R_mat[results.mask_missing] .= NaN
    end

    has_dates = !isempty(t) && isa(first(t), Dates.AbstractDateTime)
    x = has_dates ? t : (1:length(t))
    xlabel = has_dates ? "Time" : "Time Index / Step"

    p = heatmap(
        x, z, R_mat,
        clims=(-2.0, 2.0),
        color=:coolwarm,
        xlabel=xlabel,
        ylabel="Height z (m)",
        title="Dynamic Nocturnal GSPT Transition Surface (R_coord)",
        colorbar_title="R_coord",
        size=(1000, 500)
    )

    savefig(p, save_path)
    return p
end

end # module GSPTTimeLoop