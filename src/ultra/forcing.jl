using .UltraForcing

# Example metadata wrapper matching CoreTypes.ProfileMetadata
struct ProfileMetadata
    friction_velocity::Float64
    obukhov_length::Float64
    reference_height::Float64
end

"""
    build_forcing_from_smear(data::Dict{String, Float64}, nz::Int; ref_height=16.8)

Extracts SMEAR station data values and generates SCM forcing object.
"""
function build_forcing_from_smear(
    data::Dict{String,Float64},
    nz::Int;
    ref_height::Float64=16.8,
    station_prefix::String="HYY_EDDY233"
)
    # Extract variables dynamically based on station keys
    h_flux = get(data, "$(station_prefix).H", 0.0)
    le_flux = get(data, "$(station_prefix).LE", 0.0)
    u_wind = get(data, "$(station_prefix).U", 1.0)
    u_star = get(data, "$(station_prefix).u_star", 0.1)
    l_mo = get(data, "$(station_prefix).MO_length", 1e5)

    meta = ProfileMetadata(u_star, l_mo, ref_height)

    return generate_scm_forcing(
        meta,
        h_flux,
        le_flux,
        u_wind,
        nz;
        prescribed_surface_fluxes=true
    )
end