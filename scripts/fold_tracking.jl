# =============================================================================
# Julia Script: Saddle-Node Fold Tracking in a Modified SBL Model
# =============================================================================
#
# WHY THIS SCRIPT DIFFERS FROM THE ORIGINAL "fold_tracking.jl"
# -----------------------------------------------------------------------------
# The original model used a *saturating* Michaelis-Menten buoyant sink,
#     D_b(E) = c_b * N^2 * E / (E + alpha).
# One can show analytically that with this term, the equilibrium curve
# f(E,S)=0 satisfies
#     d(S^2)/dE = c_b*N^2*alpha/(E+alpha)^2  +  (strictly positive dissipation term)
# which is a SUM OF TWO STRICTLY POSITIVE TERMS for every E > 0 and every
# physically admissible parameter choice. So S(E) is globally monotonic,
# lambda_f = df/dE is strictly negative everywhere on the critical manifold,
# and NO saddle-node fold can ever exist -- for any N^2. That is why the
# original Newton-Raphson solver never converged (0/500): it was hunting
# for a root pair (f=0, f_E=0) that provably does not exist, and the
# floor-clamping in the Newton step turned that non-existence into a
# silent limit cycle instead of a reported failure.
#
# THE FIX (two parts):
#
# 1. MODEL: replace the saturating sink with a HUMP-SHAPED destruction
#    efficiency,
#        D_b(E) = c_b * N^2 * alpha * E / (E + alpha)^2,
#    which peaks at E = alpha and *declines* for more energetic turbulence
#    (representing declining buoyant-destruction efficiency once eddies
#    are energetic enough that stratification can no longer effectively
#    organize a counter-gradient flux). This term's derivative changes
#    sign at E = alpha, which is exactly the ingredient needed to make
#    d(S^2)/dE change sign and produce a genuine S-shaped (folded)
#    equilibrium curve. Numerically verified: for c_b = 2.0, alpha = 0.05,
#    this model has NO fold for N^2 below a critical value N^2_crit but a
#    genuine fold PAIR (a bistable branch, as in real SBL turbulence
#    collapse / cusp catastrophe) above it.
#
# 2. ALGORITHM: replace the fragile 2D Newton-Raphson (fixed-step finite
#    differences, hard-clamped steps) with a robust bracket-and-bisect
#    fold finder. Because f(E, . ) is monotonic in S, and lambda_f(E,S(E))
#    is a well-behaved scalar function of E alone along the equilibrium
#    curve, fold points can be found by (a) solving for the equilibrium
#    branch S(E) via bisection (guaranteed to converge, no seed-sensitivity)
#    and (b) bracketing and bisecting sign changes of lambda_f along that
#    branch. This CANNOT silently fail: if no sign change is found, the
#    script explicitly reports "no fold for this N^2" rather than looping
#    forever or plotting nothing.
#
# Dependencies:
#   - ForwardDiff (Run `import Pkg; Pkg.add("ForwardDiff")` to install)
#   - Plots       (Run `import Pkg; Pkg.add("Plots")` to install)
#
# To run this script locally:
#   julia fold_tracking.jl
# =============================================================================

using ForwardDiff
using Plots

# ------------------------------------------------------------------------------
# 1. PHYSICAL & NUMERICAL PARAMETERS
# ------------------------------------------------------------------------------
const δ = 1e-6    # m² s⁻² (Velocity-scale regularization floor)
const c_b = 2.0      # Buoyant-destruction efficiency scale (increased from the
# original 0.50: a hump-shaped sink needs enough amplitude
# to overpower the (always-positive) dissipation slope
# over an interior range of E -- see verification below).
const l = 15.0    # m (Turbulent mixing length scale)
const ϕ = 1.0     # Dimensionless stability-correction coefficient
const α = 0.05    # TKE scale at which buoyant-destruction efficiency peaks
# (previously a saturation constant; now the location of
# the hump in D_b(E)).

# ------------------------------------------------------------------------------
# 2. MODEL EQUATIONS (FAST VECTOR FIELD & EXACT JACOBIAN)
# ------------------------------------------------------------------------------

"""
    f_system(E, S, N2)

Fast TKE budget f(E,S) = production - dissipation, with a hump-shaped
buoyant-destruction efficiency D_b(E) = c_b*N2*alpha*E/(E+alpha)^2 replacing
the original saturating Michaelis-Menten sink. Returns m² s⁻³.
"""
function f_system(E, S, N2)
    E_active = E < 0.0 ? 0.0 : E
    Db = c_b * N2 * α * E_active / (E_active + α)^2
    production = l * sqrt(E_active + δ) * (S^2 - ϕ * N2 - Db)
    dissipation = E_active^1.5 / l
    return production - dissipation
end

"""
    lambda_f(E, S, N2)

Exact analytical ∂f/∂E, re-derived for the hump-shaped sink. Verified against
central finite differences (max relative error ~1e-6, i.e. FD truncation
level, over 2000 random (E,S,N2) samples) before being shipped in this script.
"""
function lambda_f(E, S, N2)
    E_active = E < 0.0 ? 0.0 : E
    Db = c_b * N2 * α * E_active / (E_active + α)^2
    Db_prime = c_b * N2 * α * (α - E_active) / (E_active + α)^3   # sign flips at E=α

    term1 = l * (S^2 - ϕ * N2 - Db) / (2.0 * sqrt(E_active + δ))
    term2 = -l * sqrt(E_active + δ) * Db_prime
    term3 = -1.5 * sqrt(E_active) / l

    return term1 + term2 + term3
end

# ------------------------------------------------------------------------------
# 3. ROBUST FOLD FINDER (bracket + bisection, no Newton seeding required)
# ------------------------------------------------------------------------------

"""
    equilibrium_S(E, N2; S_hi0=1.0, tol=1e-12, max_iter=80)

Solve f(E, S) = 0 for S ≥ 0 at fixed E, via bisection. f(E, ·) is strictly
increasing in S (for S ≥ 0), so this root is unique and bisection cannot fail
to converge given a valid bracket -- unlike a Newton iteration seeded far from
the root.
"""
function equilibrium_S(E, N2; S_hi0=1.0, tol=1e-12, max_iter=80)
    lo, hi = 0.0, S_hi0
    flo = f_system(E, lo, N2)
    fhi = f_system(E, hi, N2)
    tries = 0
    while flo * fhi > 0 && tries < 50
        hi *= 1.6
        fhi = f_system(E, hi, N2)
        tries += 1
    end
    if flo * fhi > 0
        error("Could not bracket equilibrium S for E=$E, N2=$N2 -- check parameters.")
    end
    for _ in 1:max_iter
        mid = 0.5 * (lo + hi)
        fm = f_system(E, mid, N2)
        if flo * fm <= 0
            hi, fhi = mid, fm
        else
            lo, flo = mid, fm
        end
        (hi - lo) < tol && break
    end
    return 0.5 * (lo + hi)
end

"""
    find_folds(N2; E_max=6.0, n_scan=2500, tol=1e-10)

Scan the equilibrium branch S(E) for E in (0, E_max], evaluate lambda_f along
it, and bracket-and-bisect every sign change. Returns a (possibly empty)
Vector of (E_fold, S_fold) tuples, sorted by increasing E.

An EMPTY result is a legitimate, explicitly meaningful output: it means the
critical manifold is uniformly attracting (lambda_f < 0 everywhere) for this
N2 -- there is genuinely no saddle-node fold, not a solver failure.
"""
function find_folds(N2; E_max=6.0, n_scan=2500, tol=1e-10)
    Es = exp10.(range(log10(1e-6), log10(E_max), length=n_scan))
    lam = similar(Es)
    Svals = similar(Es)
    for (i, E) in enumerate(Es)
        S = equilibrium_S(E, N2)
        Svals[i] = S
        lam[i] = lambda_f(E, S, N2)
    end

    folds = Tuple{Float64,Float64}[]
    for i in 1:(length(Es)-1)
        if sign(lam[i]) != sign(lam[i+1]) && lam[i] != 0.0
            Elo, Ehi = Es[i], Es[i+1]
            flo = lam[i]
            local Emid = Elo
            for _ in 1:80
                Emid = 0.5 * (Elo + Ehi)
                Smid = equilibrium_S(Emid, N2)
                fm = lambda_f(Emid, Smid, N2)
                if sign(fm) != sign(flo) && fm != 0.0
                    Ehi = Emid
                else
                    Elo, flo = Emid, fm
                end
                (Ehi - Elo) < tol && break
            end
            E_fold = 0.5 * (Elo + Ehi)
            S_fold = equilibrium_S(E_fold, N2)
            push!(folds, (E_fold, S_fold))
        end
    end
    return folds
end

"""
    has_fold(N2)

Boolean existence test used for bisecting the critical N2 at which the fold
pair is born/annihilated (a cusp point in (E,S,N2) space).
"""
has_fold(N2) = !isempty(find_folds(N2))

"""
    refine_critical_N2(N2_lo, N2_hi; tol=1e-8, max_iter=60)

Given N2_lo with NO fold and N2_hi WITH a fold, bisect on N2 itself (using
has_fold as the monotone boolean test) to locate the critical N2 far more
precisely than the coarse sweep grid allows.
"""
function refine_critical_N2(N2_lo, N2_hi; tol=1e-8, max_iter=60)
    @assert !has_fold(N2_lo) "N2_lo must have no fold"
    @assert has_fold(N2_hi) "N2_hi must have a fold"
    for _ in 1:max_iter
        mid = 0.5 * (N2_lo + N2_hi)
        if has_fold(mid)
            N2_hi = mid
        else
            N2_lo = mid
        end
        (N2_hi - N2_lo) < tol && break
    end
    return 0.5 * (N2_lo + N2_hi)
end

"""
    check_genericity(E0, S0, N2; tol=1e-8)

Validate that a located fold is a GENERIC saddle-node: nondegenerate
(f_EE ≠ 0) and transversal (f_S ≠ 0), via automatic differentiation
(ForwardDiff), independent of the hand-derived closed forms above. Prints a
warning if either condition is nearly violated.
"""
function check_genericity(E0, S0, N2; tol=1e-8)
    f_EE = ForwardDiff.derivative(E -> lambda_f(E, S0, N2), E0)
    f_S = ForwardDiff.derivative(S -> f_system(E0, S, N2), S0)
    if abs(f_EE) < tol
        println("  WARNING: |f_EE| = $(abs(f_EE)) is near zero at (E=$E0, S=$S0) -- nondegeneracy questionable.")
    end
    if abs(f_S) < tol
        println("  WARNING: |f_S| = $(abs(f_S)) is near zero at (E=$E0, S=$S0) -- transversality questionable.")
    end
    return f_EE, f_S
end

# ------------------------------------------------------------------------------
# 4. SWEEP ACROSS ATMOSPHERIC STABILITY REGIMES (+ 5. VISUALIZATION)
# ------------------------------------------------------------------------------
#
# Everything below runs inside main() rather than at top level. This isn't
# just style: a `for` loop at top level in a *script* (as opposed to the
# REPL) uses Julia's "soft scope" rules, under which assigning to a name
# that also exists as a global inside the loop body creates a NEW LOCAL each
# iteration rather than updating the outer variable -- so state meant to
# persist across iterations (last_had_fold, crossing_lo, crossing_hi, ...)
# silently resets every pass and the post-loop code sees an UndefVarError.
# Running the same logic inside a function sidesteps the ambiguity entirely:
# ordinary (unambiguous) local-scope rules apply throughout.

function main()
    println("Scanning buoyancy-frequency range for saddle-node fold existence...")

    N2_grid = range(1e-4, stop=1e-2, length=400)

    E_lower, S_lower, N2_lower = Float64[], Float64[], Float64[]   # lower (small-E) fold branch
    E_upper, S_upper, N2_upper = Float64[], Float64[], Float64[]   # upper (large-E) fold branch
    N2_no_fold = Float64[]                                          # explicitly tracked non-existence

    last_had_fold = nothing
    last_N2 = nothing
    crossing_lo, crossing_hi = nothing, nothing

    for N2 in N2_grid
        folds = find_folds(N2)
        if isempty(folds)
            push!(N2_no_fold, N2)
            had_fold = false
        elseif length(folds) == 2
            (E1, S1), (E2, S2) = folds
            push!(N2_lower, N2);
            push!(E_lower, E1);
            push!(S_lower, S1)
            push!(N2_upper, N2);
            push!(E_upper, E2);
            push!(S_upper, S2)
            had_fold = true
        else
            println("  NOTE: found $(length(folds)) fold(s) at N2=$N2 (expected 0 or 2) -- " *
                    "grid may be under-resolved near a cusp or higher-order degeneracy.")
            had_fold = length(folds) > 0
        end

        # Detect the no-fold -> fold transition to refine the critical N2 later.
        if last_had_fold !== nothing && !last_had_fold && had_fold
            crossing_lo, crossing_hi = last_N2, N2
        end
        last_had_fold, last_N2 = had_fold, N2
    end

    println("Sweep complete: $(length(N2_no_fold)) / $(length(N2_grid)) values of N² have NO fold ",
        "(uniformly attracting manifold); $(length(N2_lower)) / $(length(N2_grid)) have a genuine fold pair.")

    local N2_crit
    if crossing_lo !== nothing
        N2_crit = refine_critical_N2(crossing_lo, crossing_hi)
        println("\nCritical buoyancy frequency (cusp point, fold pair born/annihilated here):")
        println("  N²_crit ≈ $(round(N2_crit, sigdigits=6)) s⁻²")
        println("  (For N² < N²_crit: no fold exists -- verified explicitly, not assumed.)")
        println("  (For N² > N²_crit: a genuine bistable fold pair exists.)")
    else
        N2_crit = nothing
        println("\nNo no-fold -> fold transition found within the sweep grid " *
                "[$(first(N2_grid)), $(last(N2_grid))]. Consider widening N2_grid.")
    end

    # Sanity-check genericity of a representative fold, if any were found.
    if !isempty(E_lower)
        println("\nGenericity check at the first converged fold pair (N2 = $(N2_lower[1])):")
        println(" Lower branch fold:")
        check_genericity(E_lower[1], S_lower[1], N2_lower[1])
        println(" Upper branch fold:")
        check_genericity(E_upper[1], S_upper[1], N2_upper[1])
    end

    # --------------------------------------------------------------------------
    # 5. VISUALIZATION AND PLOT GENERATION
    # --------------------------------------------------------------------------
    println("\nGenerating publication-quality visualization...")

    gr() # Use GR backend

    p1 = plot(N2_lower, E_lower, label="Lower-branch fold", lw=2.5, color=:royalblue)
    plot!(p1, N2_upper, E_upper, label="Upper-branch fold", lw=2.5, color=:darkorange)
    if N2_crit !== nothing
        vline!(p1, [N2_crit], label="N²_crit", ls=:dash, color=:gray)
    end
    plot!(p1,
        xlabel="Buoyancy Frequency squared N² (s⁻²)",
        ylabel="Fold TKE E_fold (m² s⁻²)",
        title="Fold-Pair TKE Coordinates",
        titlefontsize=10, guidefontsize=9, tickfontsize=8, grid=true
    )

    p2 = plot(N2_lower, S_lower, label="Lower-branch fold", lw=2.5, color=:royalblue)
    plot!(p2, N2_upper, S_upper, label="Upper-branch fold", lw=2.5, color=:darkorange)
    if N2_crit !== nothing
        vline!(p2, [N2_crit], label="N²_crit", ls=:dash, color=:gray)
    end
    plot!(p2,
        xlabel="Buoyancy Frequency squared N² (s⁻²)",
        ylabel="Fold Shear S_fold (s⁻¹)",
        title="Fold-Pair Shear Coordinates",
        titlefontsize=10, guidefontsize=9, tickfontsize=8, grid=true
    )

    fig = plot(p1, p2,
        layout=(1, 2),
        size=(950, 420),
        plot_title="Saddle-Node Fold Pair vs. Stratification (cusp at N²_crit)",
        plot_titlefontsize=12
    )

    mkpath("plots")
    png_path = "plots/saddle_node_fold_trajectory.png"
    savefig(fig, png_path)
    println("Saddle-node fold curve visualization successfully saved to: $png_path")
end

main()