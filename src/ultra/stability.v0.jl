using .UltraForcing
using .UltraStability

"""
    process_smear_frame(data::Dict{String, Float64}, nz::Int; z_ref::Float64=16.8)

Evaluates forcing variables and assigns the atmospheric stability class.
"""
function process_smear_frame(
    data::Dict{String,Float64},
    nz::Int;
    z_ref::Float64=16.8,
    station_prefix::String="HYY_EDDY233"
)
    # Extract SMEAR eddy covariance variables
    u_star = get(data, "$(station_prefix).u_star", 0.1)
    l_mo = get(data, "$(station_prefix).MO_length", 1e5)
    h_flux = get(data, "$(station_prefix).H", 0.0)
    le_flux = get(data, "$(station_prefix).LE", 0.0)
    u_wind = get(data, "$(station_prefix).U", 1.0)

    # 1. Build forcing dictionary
    meta = ProfileMetadata(u_star, l_mo, z_ref)
    forcing = generate_scm_forcing(meta, h_flux, le_flux, u_wind, nz)

    # 2. Classify stability state from computed zeta
    regime = classify_stability(forcing.zeta_reference)

    return forcing, regime
end