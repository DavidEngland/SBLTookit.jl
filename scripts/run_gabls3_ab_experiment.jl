using Printf
using Statistics
using SBLToolkit

const HOURS = 24
const SECONDS_PER_HOUR = 3600

function run_experiment(closure_mode::Symbol)
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
        push!(diagnostics, step_scm!(state, config, gating_params, gs_config;
            closure_mode, enable_gating=closure_mode !== :control))
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
    override_steps = sum(count.(getfield.(diagnostics, :override_mask)))
    @printf("%s\n", label)
    @printf("  final T_2m: %.3f K\n", t2m[end])
    @printf("  mean T_2m: %.3f K\n", mean(t2m))
    @printf("  mean z_i: %.2f m\n", mean(z_i))
    @printf("  mean H_0: %.2f W m^-2\n", mean(surface_flux))
    @printf("  mean R_n: %.2f W m^-2\n", mean(radiation))
    @printf("  mean G_s: %.2f W m^-2\n", mean(ground_flux))
    @printf("  gated level-timesteps: %d\n", gated_steps)
    @printf("  eligible Ri-shutdown overrides: %d\n", override_steps)
    return t2m[end]
end

control = run_experiment(:control)
gate_only = run_experiment(:gate_only)
full_gspt = run_experiment(:full_gspt)
control_t2m = summarize("Experiment A: control", control)
gate_only_t2m = summarize("Experiment B: gate-only", gate_only)
full_gspt_t2m = summarize("Experiment C: full GSPT", full_gspt)
@printf("Final T_2m difference (gate-only - control): %.3f K\n", gate_only_t2m - control_t2m)
@printf("Final T_2m difference (full GSPT - control): %.3f K\n", full_gspt_t2m - control_t2m)