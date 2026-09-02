module SmoothingSplines

import StatsBase: fit!, fit, RegressionModel, rle, ordinalrank, mean, predict, Weights
using Reexport
using LinearAlgebra

export SmoothingSpline

@reexport using StatsBase

const LAPACKFloat = Union{Float32, Float64}

include("matrices.jl")

mutable struct SmoothingSpline{T<:LAPACKFloat} <: RegressionModel
    Xorig::Vector{T}   # original grid points, sorted
    Yorig::Vector{T}   # original values, sorted according to x
    Xrank::Vector{Int}
    Xdesign::Vector{T}
    Xcount::Vector{Int} # observation count per unique X
    Ydesign::Vector{T}
    weights::Vector{T}
    RpαQtQ::Matrix{T}  # symmetric banded matrix format
    g::Vector{T}       # fitted values
    γ::Vector{T}       # 2nd derivatives of fitted values
    λ::T
end

function fit(::Type{SmoothingSpline}, X::AbstractVector{T}, Y::AbstractVector{T}, λ::T, wts::AbstractVector{T}=fill(one(T), length(Y))) where T <: LAPACKFloat
    Xrank = ordinalrank(X)
    Xperm = sortperm(X)
    Xorig = X[Xperm]
    Yorig = Y[Xperm]

    Xdesign, Xcount = rle(Xorig)
    ws = zero(Xdesign)
    Ydesign = zero(Xdesign)
    running_rle_mean!(Ydesign, ws, Yorig, Xcount, wts[Xperm])

    RpαQtQ = QtQpR(diff(Xdesign), λ, ws)
    pbtrf!('U', 2, RpαQtQ)

    spl = SmoothingSpline{T}(
        Xorig, Yorig, Xrank, Xdesign, Xcount, Ydesign, ws, RpαQtQ,
        zero(Xdesign), fill(zero(T), length(Xdesign) - 2), λ
    )
    fit!(spl)
end

function fit!(spl::SmoothingSpline{T}) where T <: LAPACKFloat
    Y = spl.Ydesign
    ws = spl.weights
    g = spl.g
    h = diff(spl.Xdesign)
    λ = spl.λ
    Q = ReinschQ(h)

    RpαQtQ = spl.RpαQtQ
    γ = mul!(spl.γ, transpose(Q), Y)
    pbtrs!('U', 2, RpαQtQ, γ)

    mul!(g, Q, γ)
    g .= Y .- (λ .* g ./ ws)
    spl
end

function fit!(spl::SmoothingSpline{T}, Y::AbstractVector{T}) where T <: LAPACKFloat
    spl.Yorig .= Y[sortperm(spl.Xorig)]
    running_rle_mean!(spl.Ydesign, spl.weights, spl.Yorig, spl.Xcount, spl.weights)
    fit!(spl)
end

function predict(spl::SmoothingSpline{T}) where T <: LAPACKFloat
    Xcount = spl.Xcount
    curridx = 1
    g = Vector{T}(undef, length(spl.Yorig))

    @inbounds for i in eachindex(Xcount)
        for j in 1:Xcount[i]
            g[curridx] = spl.g[i]
            curridx += 1
        end
    end
    g[spl.Xrank]
end

function predict(spl::SmoothingSpline{T}, x::T) where T <: LAPACKFloat
    n = length(spl.Xdesign)
    idxl = searchsortedlast(spl.Xdesign, x)
    idxr = idxl + 1

    if idxl == 0 # Linear extrapolation left
        gl = spl.g[1]
        gr = spl.g[2]
        γ  = spl.γ[1]
        xl = spl.Xdesign[1]
        xr = spl.Xdesign[2]
        gprime = (gr - gl) / (xr - xl) - (xr - xl) * γ / 6
        return gl - (xl - x) * gprime

    elseif idxl == n # Linear extrapolation right
        gl = spl.g[n - 1]
        gr = spl.g[n]
        γ  = spl.γ[n - 2]
        xl = spl.Xdesign[n - 1]
        xr = spl.Xdesign[n]
        gprime = (gr - gl) / (xr - xl) + (xr - xl) * γ / 6
        return gr + (x - xr) * gprime

    else # Cubic interpolation
        xl = spl.Xdesign[idxl]
        xr = spl.Xdesign[idxr]
        γl = idxl == 1 ? zero(T) : spl.γ[idxl - 1]
        γr = idxl == n - 1 ? zero(T) : spl.γ[idxr - 1]
        gl = spl.g[idxl]
        gr = spl.g[idxr]
        h = xr - xl

        val = ((x - xl) * gr + (xr - x) * gl) / h
        val -= (x - xl) * (xr - x) * ((1 + (x - xl) / h) * γr + (1 + (xr - x) / h) * γl) / 6
        return val
    end
end

function predict(spl::SmoothingSpline{T}, xs::AbstractVector{T}) where T <: LAPACKFloat
    g = similar(xs)
    @inbounds for i in eachindex(xs)
        g[i] = predict(spl, xs[i])
    end
    g
end

function running_rle_mean!(g::AbstractVector{T}, w::AbstractVector{T}, Y::AbstractVector{T}, rlecount::AbstractVector{Int}, ws::AbstractVector{T}) where T <: Real
    length(g) == length(rlecount) || throw(DimensionMismatch("g and rlecount mismatch"))
    length(Y) == length(ws) || throw(DimensionMismatch("Y and ws mismatch"))

    curridx = 1
    for i in eachindex(rlecount)
        idxrange = curridx:(curridx + rlecount[i] - 1)
        g[i] = mean(Y[idxrange], Weights(ws[idxrange]))
        w[i] = sum(ws[idxrange])
        curridx += rlecount[i]
    end
    g
end

end # module