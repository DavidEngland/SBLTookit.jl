# src/GSPTTimeLoop.jl
module GSPTTimeLoop

using NCDatasets, DynamicQuantities, Plots, Dates
using ..GSPTPhase2

export ingest_netcdf_gspt, plot_gspt_transition

"""
    get_nc_var(ds::NCDataset, candidate_names::Vector{String}; required=true)

Searches the NetCDF dataset `ds` for the first matching key from `candidate_names`.
"""
function get_nc_var(ds::NCDataset, candidate_names::Vector{String}; required=true)
    for key in candidate_names
        if haskey(ds, key)
            return Array{Float64}(ds[key][:])
        end
    end
    if required
        available = join(keys(ds), ", ")
        error("None of the keys $(candidate_names) found in NetCDF dataset. Available keys: [$available]")
    end
    return nothing
end

"""
    ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
"""
function ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
    NCDataset(nc_path, "r") do ds
        # 1. Flexible Variable Resolution
        z = get_nc_var(ds, ["z", "height", "lev", "level", "alt", "heights"])
        time_raw = get_nc_var(ds, ["time", "t", "datetime"])

        N_z = length(z)
        N_t = length(time_raw)

        # 2. Field Resolution
        u_mat = get_nc_var(ds, ["u", "U", "u_wind", "wind_u", "eastward_wind"])
        v_mat = get_nc_var(ds, ["v", "V", "v_wind", "wind_v", "northward_wind"])
        th_mat = get_nc_var(ds, ["theta_v", "theta", "th", "THETA", "thv", "temp"])
        wth_mat = get_nc_var(ds, ["wthv", "wth", "w_theta", "w_th", "flux_th", "wt", "heat_flux"])
        uw_mat = get_nc_var(ds, ["uw", "u_w", "uw_flux", "momentum_flux_u"])
        vw_mat = get_nc_var(ds, ["vw", "v_w", "vw_flux", "momentum_flux_v"])

        # Optional fields
        tke_mat = get_nc_var(ds, ["tke", "TKE", "q2"]; required=false)
        Km_mat = get_nc_var(ds, ["Km", "km", "KM", "K_m"]; required=false)

        tke_mat = tke_mat === nothing ? zeros(N_z, N_t) : tke_mat
        Km_mat = Km_mat === nothing ? zeros(N_z, N_t) : Km_mat

        # Ensure correct dimension orientation (N_z, N_t)
        if size(u_mat, 1) != N_z && size(u_mat, 2) == N_z
            u_mat = permutedims(u_mat, (2, 1))
            v_mat = permutedims(v_mat, (2, 1))
            th_mat = permutedims(th_mat, (2, 1))
            wth_mat = permutedims(wth_mat, (2, 1))
            uw_mat = permutedims(uw_mat, (2, 1))
            vw_mat = permutedims(vw_mat, (2, 1))
            tke_mat = permutedims(tke_mat, (2, 1))
            Km_mat = permutedims(Km_mat, (2, 1))
        end

        # 3. Pre-allocate Output Matrices
        R_coord_2d = fill(NaN, N_z, N_t)
        C_const_2d = fill(NaN, N_z, N_t)
        C_coord_2d = fill(NaN, N_z, N_t)
        delta_clos_2d = fill(NaN, N_z, N_t)
        mask_2d = fill(false, N_z, N_t)

        # 4. Temporal Execution Loop
        for t in 1:N_t
            # Nocturnal filter: surface kinematic heat flux < 0
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
        title="Dynamic Nocturnal GSPT Transition Surface (\\mathcal{R}_{coord})",
        colorbar_title="\\mathcal{R}_{coord}",
        size=(1000, 500)
    )

    savefig(p, save_path)
    return p
end

end # module GSPTTimeLoop