module UltraStability

export classify_stability

"""
    classify_stability(zeta::Float64) -> Symbol

Categorizes boundary layer stability into regime symbols based on
the Monin-Obukhov stability parameter ζ = z/L.
"""
function classify_stability(zeta::Float64)::Symbol
    isnan(zeta) && return :unknown
    zeta > 1.0 && return :strongly_stable
    zeta > 0.1 && return :stable
    zeta >= -0.1 && return :near_neutral
    zeta >= -1.0 && return :unstable
    return :strongly_unstable
end

end # module UltraStability