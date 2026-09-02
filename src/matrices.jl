using LinearAlgebra
using LinearAlgebra.BLAS: @blasfunc
using LinearAlgebra.LAPACK: liblapack, BlasInt

# Custom safety checks replacing internal LinearAlgebra functions
function check_uplo(uplo::Char)
    (uplo == 'U' || uplo == 'L') || throw(ArgumentError("uplo argument must be 'U' or 'L', got $uplo"))
end

function check_stride1(A...)
    for a in A
        stride(a, 1) == 1 || throw(ArgumentError("array must have unit stride"))
    end
end

# --- ReinschQ Matrix ---
struct ReinschQ{T} <: AbstractMatrix{T}
    h::Vector{T}
end

function Base.size(Q::ReinschQ)
    n = length(Q.h)
    (n + 1, n - 1)
end

function Base.getindex(Q::ReinschQ{T}, i::Int, j::Int) where T
    @boundscheck checkbounds(Q, i, j)
    h = Q.h
    if i == j
        inv(h[i])
    elseif i == j + 1
        -inv(h[j]) - inv(h[j+1])
    elseif i == j + 2
        inv(h[j+1])
    else
        zero(T)
    end
end

function LinearAlgebra.mul!(out::AbstractVector, Q::Transpose{T, ReinschQ{T}}, g::AbstractVector) where T
    n = length(out)
    h = Q.parent.h
    n == size(Q, 1) || throw(DimensionMismatch("out length mismatch"))
    length(g) == size(Q, 2) || throw(DimensionMismatch("g length mismatch"))

    Δgp1 = (g[2] - g[1]) / h[1]
    @inbounds for i in 1:length(out)
        Δg = Δgp1
        Δgp1 = (g[i+2] - g[i+1]) / h[i+1]
        out[i] = Δgp1 - Δg
    end
    out
end

function LinearAlgebra.mul!(out::AbstractVector{T}, Q::ReinschQ{T}, g::AbstractVector{T}) where T
    n = length(out)
    n == size(Q, 1) || throw(DimensionMismatch("out length mismatch"))
    length(g) == size(Q, 2) || throw(DimensionMismatch("g length mismatch"))

    @inbounds for i in 1:n
        out[i] = zero(T)
        for j in max(1, i - 2):min(i, n - 2)
            out[i] += g[j] * Q[i, j]
        end
    end
    out
end

# --- ReinschR Matrix ---
struct ReinschR{T} <: AbstractMatrix{T}
    h::Vector{T}
end

function Base.size(R::ReinschR)
    n = length(R.h)
    (n - 1, n - 1)
end

function Base.getindex(R::ReinschR{T}, i::Int, j::Int) where T
    @boundscheck checkbounds(R, i, j)
    h = R.h
    if i == j
        (h[i] + h[i+1]) / 3
    elseif abs(i - j) == 1
        h[max(i, j)] / 6
    else
        zero(T)
    end
end

# --- QtQpR System Construction ---
function QtQpR(h::AbstractVector{T}, α::T, w::AbstractVector{T}=ones(T, length(h)+1)) where T <: Real
    n = length(h) - 1
    Q = ReinschQ(h)
    R = ReinschR(h)
    out = zeros(T, 3, n)

    # Main diagonal
    for i in 1:n
        out[3, i] = α * ((Q[i+2, i] / w[i+2] - Q[i+1, i] / w[i+1]) / h[i+1] -
                         (Q[i+1, i] / w[i+1] - Q[i, i] / w[i]) / h[i]) + R[i, i]
    end

    # 1st superdiagonal
    for i in 1:(n - 1)
        out[2, i+1] = α * ((Q[i+2, i+1] / w[i+2] - Q[i+1, i+1] / w[i+1]) / h[i+1] -
                           Q[i+1, i+1] / (w[i+1] * h[i])) + R[i, i+1]
    end

    # 2nd superdiagonal
    for i in 1:(n - 2)
        out[1, i+2] = α * (Q[i+2, i+2] / (w[i+2] * h[i+1])) + R[i, i+2]
    end

    out
end

# --- LAPACK Ccall Wrappers ---
for (pbtrf, pbtrs, elty) in
    ((:dpbtrf_, :dpbtrs_, :Float64),
     (:spbtrf_, :spbtrs_, :Float32),
     (:zpbtrf_, :zpbtrs_, :ComplexF64),
     (:cpbtrf_, :cpbtrs_, :ComplexF32))
    @eval begin
        function pbtrf!(uplo::Char, kd::Integer, AB::StridedMatrix{$elty})
            check_uplo(uplo)
            check_stride1(AB)
            n = size(AB, 2)
            info = Ref{BlasInt}()
            ccall((@blasfunc($pbtrf), liblapack), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ref{BlasInt}),
                  UInt8(uplo), n, kd, AB, max(1, stride(AB, 2)), info)
            info[] == 0 || error("LAPACK pbtrf! failed with info = $(info[])")
            AB
        end

        function pbtrs!(uplo::Char, kd::Integer,
                        AB::StridedMatrix{$elty}, B::StridedVecOrMat{$elty})
            check_uplo(uplo)
            check_stride1(AB, B)
            info = Ref{BlasInt}()
            n = size(AB, 2)
            if n != size(B, 1)
                throw(DimensionMismatch("Matrix AB has dimensions $(size(AB)), but RHS B has dimensions $(size(B))"))
            end
            ccall((@blasfunc($pbtrs), liblapack), Cvoid,
                  (Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
                   Ptr{$elty}, Ref{BlasInt}, Ptr{$elty}, Ref{BlasInt},
                   Ref{BlasInt}),
                  UInt8(uplo), n, kd, size(B, 2), AB, max(1, stride(AB, 2)),
                  B, max(1, stride(B, 2)), info)
            info[] == 0 || error("LAPACK pbtrs! failed with info = $(info[])")
            B
        end
    end
end