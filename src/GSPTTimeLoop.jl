module GSPTTimeLoop

using NCDatasets, Plots, Dates, Statistics
using ..GSPTPhase2

export ingest_netcdf_gspt, plot_gspt_transition

# Atmospheric thermodynamic constants
const RHO_AIR = 1.225     # Reference air density [kg/m^3]
const CP_AIR = 1004.67   # Specific heat of dry air [J/(kg K)]
const RHO_CP = RHO_AIR * CP_AIR  # ~1230.7 J/(m^3 K)
const LAPSE_DRY = 0.0098  # Dry adiabatic lapse rate [K/m]

# Explicit matching of level-suffixed tower variables: variable + height + optional 'm'
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
        elseif N_t > 0 && length(vec_arr) % N_t == 0
            n_z_guess = length(vec_arr) ÷ N_t
            reshaped = reshape(vec_arr, n_z_guess, N_t)
            if n_z_guess == N_z
                mat .= reshaped
            elseif n_z_guess < N_z
                mat[1:n_z_guess, :] .= reshaped
                @warn "ensure_2d_matrix: reshaped 1D '$name' to ($n_z_guess,$N_t) and padded to N_z=$N_z with NaN."
            else
                mat .= reshaped[1:N_z, :]
                @warn "ensure_2d_matrix: reshaped 1D '$name' to ($n_z_guess,$N_t) and truncated to N_z=$N_z."
            end
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

"""
    interpolate_column_between_valid_levels!(col, z)

Linearly interpolates NaN gaps only between valid observed levels.
Values below the lowest valid level or above the highest valid level remain unchanged.
"""
function interpolate_column_between_valid_levels!(col::AbstractVector{Float64}, z::Vector{Float64})
    valid_idx = findall(i -> isfinite(col[i]) && isfinite(z[i]), eachindex(col))
    length(valid_idx) < 2 && return col

    for k in 1:(length(valid_idx) - 1)
        i1 = valid_idx[k]
        i2 = valid_idx[k + 1]
        i2 <= i1 + 1 && continue

        z1 = z[i1]
        z2 = z[i2]
        if !isfinite(z1) || !isfinite(z2) || z2 == z1
            continue
        end

        v1 = col[i1]
        v2 = col[i2]
        for i in (i1 + 1):(i2 - 1)
            if isnan(col[i])
                w = (z[i] - z1) / (z2 - z1)
                col[i] = v1 + w * (v2 - v1)
            end
        end
    end
    return col
end

"""
    interpolate_state_profiles!(mat, z)

Applies interior-only linear interpolation column-wise.
"""
function interpolate_state_profiles!(mat::Matrix{Float64}, z::Vector{Float64})
    N_t = size(mat, 2)
    for t in 1:N_t
        interpolate_column_between_valid_levels!(view(mat, :, t), z)
    end
    return mat
end

"""
    synthesize_flux_decay_profiles!(flux_mat, z; h_sbl=200.0)

Fills missing flux profile values with F(z) = F0 * max(0, 1 - z / h_sbl),
where F0 is the lowest-height finite flux in each time column.
"""
function synthesize_flux_decay_profiles!(flux_mat::Matrix{Float64}, z::Vector{Float64}; h_sbl::Float64=200.0)
    N_z, N_t = size(flux_mat)
    for t in 1:N_t
        finite_idx = findall(i -> isfinite(flux_mat[i, t]), 1:N_z)
        isempty(finite_idx) && continue

        i0 = argmin(z[finite_idx])
        F0 = flux_mat[finite_idx[i0], t]
        if !isfinite(F0)
            continue
        end

        for i in 1:N_z
            if isnan(flux_mat[i, t])
                ramp = max(0.0, 1.0 - z[i] / h_sbl)
                flux_mat[i, t] = F0 * ramp
            end
        end
    end
    return flux_mat
end

"""
    compute_triple_point_dispersion(z, Km, e; tke_fraction=0.05)

Computes non-dimensional triple-point dispersion for a single time-column profile:
- z_K: diffusivity cutoff marker
- z_e: TKE level-set marker
- z_ez: strongest TKE gradient marker
"""
function compute_triple_point_dispersion(z::Vector{Float64}, Km::Vector{Float64}, e::Vector{Float64}; tke_fraction::Float64=0.05)
    valid_mask = isfinite.(z) .& isfinite.(Km) .& isfinite.(e)
    if count(valid_mask) < 3
        return NaN
    end

    z_v = z[valid_mask]
    Km_v = Km[valid_mask]
    e_v = e[valid_mask]

    H_sbl = maximum(z_v) - minimum(z_v)
    H_sbl > 0.0 || return NaN

    e_min = minimum(e_v) + tke_fraction * (maximum(e_v) - minimum(e_v))
    idx_e = argmin(abs.(e_v .- e_min))
    z_e = z_v[idx_e]

    Km_min = minimum(Km_v) + tke_fraction * (maximum(Km_v) - minimum(Km_v))
    idx_K = argmin(abs.(Km_v .- Km_min))
    z_K = z_v[idx_K]

    de_dz = abs.(diff(e_v) ./ diff(z_v))
    z_ez_mid = (z_v[1:end-1] .+ z_v[2:end]) ./ 2.0
    z_ez = z_ez_mid[argmax(de_dz)]

    return (max(z_K, z_e, z_ez) - min(z_K, z_e, z_ez)) / H_sbl
end

# ---------------------------------------------------------------------------
# Layer 2: Physical & Thermodynamic Normalization
# ---------------------------------------------------------------------------

"""
    normalize_thermodynamics(ds::NCDataset, z::Vector{Float64}, N_t::Int)

Extracts and converts thermodynamic variables to Virtual Potential Temperature (\\theta_v).
"""
function normalize_thermodynamics(ds::NCDataset, z::Vector{Float64}, N_t::Int)
    N_z = length(z)

    # Direct Virtual Potential Temperature (\theta_v)
    thv_raw = get_nc_var(ds, ["theta_v", "thv", "THETA_V"]; required=false)
    if thv_raw !== nothing
        thv_mat = ensure_2d_matrix(thv_raw, N_z, N_t; name="theta_v")
        valid = filter(!isnan, thv_mat)
        if !isempty(valid) && mean(valid) < 100.0
            thv_mat .+= 273.15
        end
        return thv_mat
    end

    # Potential Temperature (\theta)
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

    # Absolute Temperature (T) conversion to \theta_v
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
    normalize_heat_flux(wth_mat; flux_convention=:auto, rho_cp=RHO_CP, cooling_flux_sign=:negative)

Normalizes heat flux to kinematic units (K m/s). Converts Sensible Heat Flux (W/m^2) using \\rho c_p.

`cooling_flux_sign` controls sign convention for cooling flux:
- `:negative` -> cooling is already negative (default)
- `:positive` -> cooling is positive in source data, so sign is flipped
- `:auto` -> no sign flip (audit output should be used to verify)
"""
function normalize_heat_flux(wth_mat::Matrix{Float64}; flux_convention::Symbol=:auto,
    rho_cp::Float64=RHO_CP, cooling_flux_sign::Symbol=:negative)
    wth_norm = copy(wth_mat)
    valid_vals = filter(!isnan, wth_norm)
    isempty(valid_vals) && return wth_norm

    is_sensible = flux_convention == :sensible || (flux_convention == :auto && maximum(abs.(valid_vals)) > 5.0)
    if is_sensible
        wth_norm ./= rho_cp
    end

    if cooling_flux_sign == :positive
        wth_norm .*= -1.0
    elseif cooling_flux_sign != :negative && cooling_flux_sign != :auto
        @warn "Unknown cooling_flux_sign=$(cooling_flux_sign). Expected :negative, :positive, or :auto. Leaving sign unchanged."
    end

    return wth_norm
end

"""
    normalize_momentum_flux(ds::NCDataset, N_z::Int, N_t::Int)

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
# Layer 3: Tower & Gridded Ingestion Workflows
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
                        flux_convention=:auto, rho_cp=RHO_CP, cooling_flux_sign=:negative,
                        nocturnal_only=true, use_threads=false, debug_audit=false,
                        mask_ill_conditioned_in_solver=true, synthesize_missing_fluxes=false,
                        h_sbl_diagnostic=200.0)

Ingests NetCDF observations, normalizes thermodynamics and physical fluxes, performs profile validation,
and executes GSPT over valid time columns. Returns explicit diagnostic masks.
"""
function ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05,
    flux_convention::Symbol=:auto, rho_cp::Float64=RHO_CP, cooling_flux_sign::Symbol=:negative,
    nocturnal_only::Bool=true, use_threads::Bool=false, debug_audit::Bool=false,
    mask_ill_conditioned_in_solver::Bool=true, synthesize_missing_fluxes::Bool=false,
    h_sbl_diagnostic::Float64=200.0, dispersion_qc_max::Float64=0.10)

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

            # Safely extract 1D height profile from 1D, 2D, or N-D coordinate arrays
            if ndims(z_raw) == 2
                s1, s2 = size(z_raw)
                if s2 == N_t
                    z = Float64.(z_raw[:, 1])
                elseif s1 == N_t
                    z = Float64.(z_raw[1, :])
                else
                    z = Float64.(z_raw[:, 1])
                end
            elseif ndims(z_raw) > 2
                z = Float64.(vec(z_raw[:, 1, 1]))
            else
                z_vec = Float64.(vec(z_raw))
                if length(z_vec) > N_t && N_t > 0 && length(z_vec) % N_t == 0
                    n_z_guess = length(z_vec) ÷ N_t
                    z = reshape(z_vec, n_z_guess, N_t)[:, 1]
                else
                    z = z_vec
                end
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

        # Normalize heat flux to kinematic units (K m/s) and harmonize sign convention for cooling
        wth_mat = normalize_heat_flux(wth_mat;
            flux_convention=flux_convention,
            rho_cp=rho_cp,
            cooling_flux_sign=cooling_flux_sign
        )

        if synthesize_missing_fluxes
            # Diagnostic-only reconstruction: fill state interior gaps without top extrapolation.
            interpolate_state_profiles!(u_mat, z)
            interpolate_state_profiles!(v_mat, z)
            interpolate_state_profiles!(thv_mat, z)

            # Diagnostic-only turbulence fallback: linear SBL decay from the lowest finite flux value.
            synthesize_flux_decay_profiles!(wth_mat, z; h_sbl=h_sbl_diagnostic)
            synthesize_flux_decay_profiles!(uw_mat, z; h_sbl=h_sbl_diagnostic)
            synthesize_flux_decay_profiles!(vw_mat, z; h_sbl=h_sbl_diagnostic)
        end

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

        mask_missing = fill(false, N_z, N_t)
        mask_nocturnal = fill(false, N_z, N_t)
        mask_ill_conditioned = fill(false, N_z, N_t)
        mask_dispersion_qc = fill(false, N_z, N_t)
        delta_tp = fill(NaN, N_t)

        skipped_nocturnal_filter = Threads.Atomic{Int}(0)
        skipped_missing_profiles = Threads.Atomic{Int}(0)
        skipped_dispersion_qc = Threads.Atomic{Int}(0)
        solved_columns = Threads.Atomic{Int}(0)
        solved_columns_with_ill = Threads.Atomic{Int}(0)
        masked_ill_points = Threads.Atomic{Int}(0)
        solved_partial_profile_columns = Threads.Atomic{Int}(0)

        # Minimum level guard: finite differencing in GSPTPhase2 requires N_z >= 3
        if N_z < 3
            @warn "Dataset '$nc_path' has only $N_z vertical levels ($z m). GSPT requires N_z >= 3 for spatial derivatives. Skipping solver execution."
            return (
                z=z, time=time_raw, R_coord=R_coord_2d, C_const=C_const_2d,
                C_coord=C_coord_2d, delta_closure=delta_clos_2d,
                mask_missing=fill(true, N_z, N_t), mask_nocturnal=mask_nocturnal,
                mask_ill_conditioned=mask_ill_conditioned,
                mask_dispersion_qc=mask_dispersion_qc,
                delta_tp=delta_tp
            )
        end

        compute_column! = t -> begin
            sfc_flux = wth_mat[1, t]

            is_noct = !isnan(sfc_flux) && sfc_flux < 0.0
            mask_nocturnal[:, t] .= is_noct
            if nocturnal_only && !is_noct
                Threads.atomic_add!(skipped_nocturnal_filter, 1)
                return nothing
            end

            Km_col = view(Km_mat, :, t)
            e_col = view(tke_mat, :, t)
            δ_TP = compute_triple_point_dispersion(z, Km_col, e_col)
            delta_tp[t] = δ_TP
            if !isfinite(δ_TP) || δ_TP >= dispersion_qc_max
                mask_dispersion_qc[:, t] .= true
                mask_ill_conditioned[:, t] .= true
                R_coord_2d[:, t] .= NaN
                Threads.atomic_add!(skipped_dispersion_qc, 1)
                return nothing
            end

            u_col, v_col, thv_col = u_mat[:, t], v_mat[:, t], thv_mat[:, t]
            has_missing = any(isnan, u_col) || any(isnan, v_col) || any(isnan, thv_col)
            if has_missing && !synthesize_missing_fluxes
                mask_missing[:, t] .= true
                Threads.atomic_add!(skipped_missing_profiles, 1)
                return nothing
            end

            if has_missing && synthesize_missing_fluxes
                valid_idx = findall(i -> isfinite(u_col[i]) && isfinite(v_col[i]) && isfinite(thv_col[i]) &&
                    isfinite(wth_mat[i, t]) && isfinite(uw_mat[i, t]) && isfinite(vw_mat[i, t]), 1:N_z)

                if length(valid_idx) < 3
                    mask_missing[:, t] .= true
                    Threads.atomic_add!(skipped_missing_profiles, 1)
                    return nothing
                end

                data_t = ProfileData(
                    z[valid_idx], u_col[valid_idx], v_col[valid_idx], thv_col[valid_idx],
                    wth_mat[valid_idx, t], uw_mat[valid_idx, t], vw_mat[valid_idx, t],
                    tke_mat[valid_idx, t], Km_mat[valid_idx, t],
                    σ_u, σ_v, σ_th, 0.01
                )

                res_t = compute_gspt(
                    data_t;
                    is_observation=true,
                    S2_min=S2_min,
                    mask_ill_conditioned=mask_ill_conditioned_in_solver
                )

                mask_missing[:, t] .= true
                mask_missing[valid_idx, t] .= false

                R_coord_2d[valid_idx, t] .= res_t.const_geom.R_coord
                C_const_2d[valid_idx, t] .= res_t.const_geom.C_const
                C_coord_2d[valid_idx, t] .= res_t.const_geom.C_coord
                delta_clos_2d[valid_idx, t] .= res_t.obs_diag.closure_residual
                mask_ill_conditioned[valid_idx, t] .= res_t.obs_diag.ill_conditioned_mask
                Threads.atomic_add!(solved_partial_profile_columns, 1)
            else
                data_t = ProfileData(
                    z, u_col, v_col, thv_col,
                    wth_mat[:, t], uw_mat[:, t], vw_mat[:, t],
                    tke_mat[:, t], Km_mat[:, t],
                    σ_u, σ_v, σ_th, 0.01
                )

                res_t = compute_gspt(
                    data_t;
                    is_observation=true,
                    S2_min=S2_min,
                    mask_ill_conditioned=mask_ill_conditioned_in_solver
                )

                R_coord_2d[:, t] .= res_t.const_geom.R_coord
                C_const_2d[:, t] .= res_t.const_geom.C_const
                C_coord_2d[:, t] .= res_t.const_geom.C_coord
                delta_clos_2d[:, t] .= res_t.obs_diag.closure_residual
                mask_ill_conditioned[:, t] .= res_t.obs_diag.ill_conditioned_mask
            end

            Threads.atomic_add!(solved_columns, 1)
            ill_points_t = count(res_t.obs_diag.ill_conditioned_mask)
            Threads.atomic_add!(masked_ill_points, ill_points_t)
            if ill_points_t > 0
                Threads.atomic_add!(solved_columns_with_ill, 1)
            end
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

        if debug_audit
            sfc = wth_mat[1, :]
            n_sfc_valid = count(!isnan, sfc)
            n_sfc_neg = count(x -> !isnan(x) && x < 0.0, sfc)
            n_sfc_pos = count(x -> !isnan(x) && x > 0.0, sfc)
            n_finite_r = count(isfinite, R_coord_2d)
            @info "GSPT ingest audit: $(basename(nc_path))"
                N_t=N_t
                N_z=N_z
                nocturnal_only=nocturnal_only
                flux_convention=flux_convention
                cooling_flux_sign=cooling_flux_sign
                synthesize_missing_fluxes=synthesize_missing_fluxes
                h_sbl_diagnostic=h_sbl_diagnostic
                dispersion_qc_max=dispersion_qc_max
                sfc_flux_valid=n_sfc_valid
                sfc_flux_negative=n_sfc_neg
                sfc_flux_positive=n_sfc_pos
                skipped_nocturnal_filter=skipped_nocturnal_filter[]
                skipped_missing_profiles=skipped_missing_profiles[]
                skipped_dispersion_qc=skipped_dispersion_qc[]
                solved_columns=solved_columns[]
                solved_partial_profile_columns=solved_partial_profile_columns[]
                solved_columns_with_ill_conditioning=solved_columns_with_ill[]
                masked_ill_conditioned_points=masked_ill_points[]
                finite_R_coord_points=n_finite_r

            if nocturnal_only && n_sfc_neg == 0 && n_sfc_pos > 0
                @warn "No negative surface heat-flux values after normalization; if this dataset uses positive cooling flux, rerun with cooling_flux_sign=:positive."
            end
        end

        return (
            z=z, time=time_raw, R_coord=R_coord_2d, C_const=C_const_2d,
            C_coord=C_coord_2d, delta_closure=delta_clos_2d,
            mask_missing=mask_missing, mask_nocturnal=mask_nocturnal,
            mask_ill_conditioned=mask_ill_conditioned,
            mask_dispersion_qc=mask_dispersion_qc,
            delta_tp=delta_tp
        )
    end
end

# ---------------------------------------------------------------------------
# Layer 4: Transition Diagnostics & Visualization
# ---------------------------------------------------------------------------

"""
    plot_gspt_transition(results; save_path="gspt_transition_heatmap.png",
                          mask_ill_conditioned=true, mask_missing=true,
                          plot_title="Dynamic Nocturnal GSPT Transition Surface (R_coord)")

Plots dynamic R_coord heatmap over height and time. Masks ill-conditioned or missing observations.
"""
function plot_gspt_transition(results; save_path="gspt_transition_heatmap.png",
    mask_ill_conditioned::Bool=true, mask_missing::Bool=true,
    plot_title::String="Dynamic Nocturnal GSPT Transition Surface (R_coord)")

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

    if count(isfinite, R_mat) == 0
        @warn "Heatmap contains no finite R_coord values after masking. Writing annotated placeholder plot to $(save_path)."
        p = plot(
            xlabel=xlabel,
            ylabel="Height z (m)",
            title=plot_title,
            legend=false,
            size=(1000, 500)
        )
        if !isempty(x) && !isempty(z)
            annotate!(p, x[cld(length(x), 2)], z[cld(length(z), 2)], text("No finite R_coord values", 12, :red))
        end
        savefig(p, save_path)
        return p
    end

    p = heatmap(
        x, z, R_mat,
        clims=(-2.0, 2.0),
        color=:coolwarm,
        xlabel=xlabel,
        ylabel="Height z (m)",
        title=plot_title,
        colorbar_title="R_coord",
        size=(1000, 500)
    )

    savefig(p, save_path)
    return p
end

end # module GSPTTimeLoop