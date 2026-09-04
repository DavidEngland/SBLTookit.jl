using Printf
using Statistics
using SBLToolkit

const HOURS = 24
const SECONDS_PER_HOUR = 3600

function run_experiment(; enable_gating::Bool)
    z = collect(5.0:5.0:200.0)
    state = initialize_cabauw_state(z)
    config = SCMConfig(z; dt=60.0, surface_mode=:slab, radiation_max=400.0,
        longwave_loss=80.0, radiation_peak_time=43_200.0, soil_temperature=279.0,
        ground_conductance=2.0)
    gating_params = BifurcationGatingParams(
        epsilon_on=-0.10, epsilon_off=-0.05, zeta_z_tol=1e-2, km_floor=0.1, kh_floor=0.01,
    )
    gs_config = GSPTModelConfig()
    diagnostics = SCMDiagnostics{Float64}[]

    for _ in 1:(HOURS * SECONDS_PER_HOUR ÷ round(Int, config.dt))
        push!(diagnostics, step_scm!(state, config, gating_params, gs_config; enable_gating))
    end
    return diagnostics
end

function summarize(label::String, diagnostics::Vector{SCMDiagnostics{Float64}})
    t2m = getfield.(diagnostics, :t2m)
    z_i = getfield.(diagnostics, :boundary_layer_height)
    surface_flux = getfield.(diagnostics, :surface_sensible_heat_flux)
    radiation = getfield.(diagnostics, :net_radiation)
    ground_flux = getfield.(diagnostics, :ground_heat_flux)
    gated_steps = sum(getfield.(diagnostics, :gated_levels))
    @printf("%s\n", label)
    @printf("  final T_2m: %.3f K\n", t2m[end])
    @printf("  mean T_2m: %.3f K\n", mean(t2m))
    @printf("  mean z_i: %.2f m\n", mean(z_i))
    @printf("  mean H_0: %.2f W m^-2\n", mean(surface_flux))
    @printf("  mean R_n: %.2f W m^-2\n", mean(radiation))
    @printf("  mean G_s: %.2f W m^-2\n", mean(ground_flux))
    @printf("  gated level-timesteps: %d\n", gated_steps)
    return t2m[end]
end

baseline = run_experiment(enable_gating=false)
gated = run_experiment(enable_gating=true)
baseline_t2m = summarize("Experiment A: baseline", baseline)
gated_t2m = summarize("Experiment B: geometry-aware gate", gated)
@printf("Final T_2m difference (B - A): %.3f K\n", gated_t2m - baseline_t2m)