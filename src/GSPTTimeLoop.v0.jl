# src/GSPTTimeLoop.jl
module GSPTTimeLoop

using NCDatasets, DynamicQuantities, Plots, Dates
using SBLToolkit.GSPTPhase2
export ingest_netcdf_gspt, plot_gspt_transition
"""
    ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)

Ingests a multi-level boundary layer NetCDF file, extracts time-series matrices,
runs primitive-space MDP for each time step, and returns 2D transition arrays.
"""
function ingest_netcdf_gspt(nc_path::String; S2_min=1e-3, σ_u=0.05, σ_v=0.05, σ_th=0.05)
    NCDataset(nc_path, "r") do ds
        # 1. Load Axes & Coordinates
        z = Array{Float64}(ds["height"][:])        # Height vector [N_z]
        time_raw = ds["time"][:]                    # Time vector [N_t]

        N_z = length(z)
        N_t = length(time_raw)

        # 2. Extract 2D Fields: Dimensions assumed (N_z, N_t)
        u_mat    = Array{Float64}(ds["u"][:])
        v_mat    = Array{Float64}(ds["v"][:])
        th_mat   = Array{Float64}(ds["theta_v"][:])
        wth_mat  = Array{Float64}(ds["wthv"][:])    # Kinematic heat flux
        uw_mat   = Array{Float64}(ds["uw"][:])      # u-momentum flux
        vw_mat   = Array{Float64}(ds["vw"][:])      # v-momentum flux
        tke_mat  = haskey(ds, "tke") ? Array{Float64}(ds["tke"][:]) : zeros(N_z, N_t)
        Km_mat   = haskey(ds, "Km")  ? Array{Float64}(ds["Km"][:])  : zeros(N_z, N_t)

        # 3. Pre-allocate 2D Output Surfaces
        R_coord_2d    = fill(NaN, N_z, N_t)
        C_const_2d    = fill(NaN, N_z, N_t)
        C_coord_2d    = fill(NaN, N_z, N_t)
        delta_clos_2d = fill(NaN, N_z, N_t)
        mask_2d       = fill(false, N_z, N_t)

        # 4. Temporal Execution Loop
        for t in 1:N_t
            # Surface heat flux nocturnal check: w'θv' < 0 at surface
            is_nocturnal = wth_mat[1, t] < 0.0
            !is_nocturnal && continue # Skip daytime convective profiles

            # Assemble slice profile
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

            # Evaluate MDP & GSPT profile
            res_t = compute_gspt(data_t; is_observation=true, S2_min=S2_min)

            # Store results into 2D transition matrices
            R_coord_2d[:, t]    .= res_t.const_geom.R_coord
            C_const_2d[:, t]    .= res_t.const_geom.C_const
            C_coord_2d[:, t]    .= res_t.const_geom.C_coord
            delta_clos_2d[:, t] .= res_t.obs_diag.closure_residual
            mask_2d[:, t]       .= res_t.obs_diag.ill_conditioned_mask
        end

        return (z=z, time=time_raw, R_coord=R_coord_2d, C_const=C_const_2d,
                C_coord=C_coord_2d, delta_closure=delta_clos_2d, mask=mask_2d)
    end
end

end # module GSPTTimeLoop