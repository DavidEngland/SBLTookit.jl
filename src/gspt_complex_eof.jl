#!/usr/bin/env julia
# ==============================================================================
# GSPT COMPLEX EOF (CEOF) & PHASE-GRADIENT SEPARATION PIPELINE
# Developed for Wave-Jet Separation in the Nocturnal Stable Boundary Layer (SBL)
# ==============================================================================
# This script implements:
# 1. Self-contained, dependency-free Complex Analytic Signal generation 
#    via a robust matrix Fourier Hilbert transform.
# 2. Singular Value Decomposition (SVD) of the complex signal matrix.
# 3. Phase extraction and a robust 1D Phase Unwrapping Operator.
# 4. Asymmetric, non-uniform spatial derivative operators (D1, D2) 
#    to evaluate the topological coordinate-invariant phase gradient (∂θ/∂z).
# 5. Dual-mode physical classification:
#    - Standing low-level jet structures (uniform phase, ∂θ/∂z ≈ 0)
#    - Propagating gravity waves (sloping phase gradient, ∂θ/∂z ≠ 0)
# ==============================================================================

using LinearAlgebra
using Printf

"""
    stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)

Computes the finite-difference stencil weights at a target coordinate `z0` using 
an arbitrary grid set `z_stencil` for a derivative of order `m`.
Uses a Vandermonde solver derived from local Taylor expansion.
"""
function stencil_weights(z_stencil::Vector{Float64}, z0::Float64, m::Int)
    p = length(z_stencil)
    # Construct Vandermonde system based on Taylor terms (z - z0)^(k-1) / (k-1)!
    A = [Float64(z - z0)^(k - 1) / factorial(k - 1) for k in 1:p, z in z_stencil]
    b = zeros(p)
    b[m + 1] = 1.0  # Select target derivative order
    return A \ b
end

"""
    build_operators(z::Vector{Float64})

Generates non-uniform derivative matrices D1 (first derivative) and D2 (second derivative)
on a vertical coordinate vector `z`. Boundary points utilize second-order accurate asymmetric 
one-sided stencils to ensure uniform error convergence across all levels.
"""
function build_operators(z::Vector{Float64})
    n = length(z)
    D1 = zeros(n, n)
    D2 = zeros(n, n)
    
    for i in 1:n
        if i == 1
            idx = [1, 2, 3]
        elseif i == n
            idx = [n-2, n-1, n]
        else
            idx = [i-1, i, i+1]
        end
        D1[i, idx] = stencil_weights(z[idx], z[i], 1)
        D2[i, idx] = stencil_weights(z[idx], z[i], 2)
    end
    return D1, D2
end

"""
    unwrap_phase(θ::Vector{Float64})

Performs 1D phase unwrapping along a spatial or temporal vector to prevent 
artificial 2π jumps from corrupting numerical gradients.
"""
function unwrap_phase(θ::Vector{Float64})
    n = length(θ)
    unwrapped = copy(θ)
    for i in 2:n
        diff = θ[i] - θ[i-1]
        # Map difference to [-π, π] interval
        diff = mod(diff + π, 2π) - π
        unwrapped[i] = unwrapped[i-1] + diff
    end
    return unwrapped
end

"""
    compute_hilbert_matrix(N::Int)

Constructs a self-contained, dependency-free Hilbert transform operator matrix.
Applying this matrix to a real vector of length N yields its complex analytical signal:
z(t) = x(t) + i * H[x(t)]
"""
function compute_hilbert_matrix(N::Int)
    # Construct the Discrete Fourier Transform (DFT) matrix F
    # F_jk = exp(-i * 2π * j * k / N)
    F = [exp(-im * 2π * (j-1) * (k-1) / N) for j in 1:N, k in 1:N]
    
    # Construct the frequency-domain analytic filter S
    S = zeros(N)
    if N % 2 == 0
        S[1] = 1.0
        S[div(N, 2) + 1] = 1.0
        S[2:div(N, 2)] .= 2.0
    else
        S[1] = 1.0
        S[2:div(N + 1, 2)] .= 2.0
    end
    
    # Inverse DFT matrix is the conjugate transpose scaled by 1/N
    Finv = conj(F) ./ N
    
    # The analytical signal operator H_op = F_inv * diag(S) * F
    H_op = Finv * (Diagonal(S) * F)
    return H_op
end

"""
    run_gspt_complex_eof(z::Vector{Float64}, U::Matrix{Float64})

Performs Complex Empirical Orthogonal Function (CEOF) analysis on a velocity matrix U (size N_z x N_t)
and maps Jet standing modes versus propagating Wave modes using the vertical phase-gradient (∂θ/∂z).
"""
function run_gspt_complex_eof(z::Vector{Float64}, U::Matrix{Float64})
    N_z, N_t = size(U)
    
    # 1. Compute temporal anomalies (subtract the time-mean velocity profile)
    U_mean = mean(U, dims=2)
    U_anomaly = U .- U_mean
    
    # 2. Generate Complex Analytic Signal row-wise using the self-contained Hilbert operator
    println("Constructing Hilbert transform operator of dimension (", N_t, "x", N_t, ")...")
    H_operator = compute_hilbert_matrix(N_t)
    U_complex = zeros(Complex{Float64}, N_z, N_t)
    for i in 1:N_z
        U_complex[i, :] = H_operator * U_anomaly[i, :]
    end
    
    # 3. Compute Complex Singular Value Decomposition (SVD)
    println("Computing Complex SVD...")
    svd_decomp = svd(U_complex)
    
    # Extract singular values and calculate variance portion
    S = svd_decomp.S
    variance_explained = (S .^ 2) ./ sum(S .^ 2)
    
    # 4. Build non-uniform spatial operators for the CASES-99 tower coordinates
    D1, D2 = build_operators(z)
    
    # We analyze the first two dominant modes (typically standing jet and propagating wave)
    results = Dict{Int, Any}()
    for k in 1:min(2, N_z)
        # Left singular vector corresponds to the spatial mode
        V_k = svd_decomp.U[:, k]
        
        # Calculate amplitude and raw phase
        amp = abs.(V_k)
        phase_raw = angle.(V_k)
        
        # Unwrap phase along the vertical z coordinate
        phase_unwrapped = unwrap_phase(phase_raw)
        
        # Compute topological vertical phase gradient (∂θ/∂z)
        phase_grad = D1 * phase_unwrapped
        
        # Store results
        results[k] = (
            amplitude = amp,
            phase_raw = phase_raw,
            phase_unwrapped = phase_unwrapped,
            phase_gradient = phase_grad,
            variance = variance_explained[k]
        )
    end
    
    return results
end

# Simple helper for mean calculation to stay standard-library free
mean(x; dims=1) = sum(x, dims=dims) ./ size(x, dims...)

# ==============================================================================
# DEMONSTRATION & VERIFICATION LAYER
# ==============================================================================

function main()
    # CASES-99 non-uniform tower coordinates
    z_tower = [1.5, 5.0, 10.0, 20.0, 30.0, 45.0, 55.0]
    N_z = length(z_tower)
    N_t = 144  # 24-hour cycle of 10-minute records
    time = range(0.0, 24.0, length=N_t)
    
    println("="^110)
    println("                          GSPT WAVE-JET COMPLEX SVD DECOMPOSITION DEMONSTRATION")
    println("="^110)
    
    # 1. Synthesize physical standing Jet Mode peaking near z ≈ 30m (no spatial phase propagation)
    U_jet_spatial = [10.0 * exp(-((h - 30.0) / 15.0)^2) for h in z_tower]
    U_jet_temporal = [1.0 + sin(2 * π * t / 24.0 - π/2) for t in time]
    U_jet = [U_jet_spatial[i] * U_jet_temporal[t] for i in 1:N_z, t in 1:N_t]
    
    # 2. Synthesize propagating Internal Gravity Wave (IGW) with vertical wavenumber m_wave
    m_wave = 0.15          # rad/m (wave propagates vertically)
    omega_wave = 2 * π / 2.0 # 2-hour oscillation period
    U_wave = [1.5 * cos(omega_wave * time[t] - m_wave * z_tower[i]) for i in 1:N_z, t in 1:N_t]
    
    # 3. Add random sonic anemometer noise (σ = 0.05 m/s)
    noise = 0.05 .* randn(N_z, N_t)
    
    # Combine fields
    U_total = U_jet .+ U_wave .+ noise
    
    # Run the CEOF separation pipeline
    results = run_gspt_complex_eof(z_tower, U_total)
    
    # Display the results
    println("\n" * "SVD MODE CLASSIFICATION STATISTICS:")
    println("-"^110)
    for k in 1:2
        @printf("Mode %d Explained Variance: %6.2f%%\n", k, results[k].variance * 100.0)
    end
    println("-"^110)
    
    for k in 1:2
        println("\n" * "="^110)
        @printf("                                      CEOF MODE %d PROFILE ANALYSIS\n", k)
        println("="^110)
        @printf("%-6s | %-15s | %-20s | %-20s | %-20s\n",
                "z (m)", "Amplitude", "Raw Phase (rad)", "Unwrapped Phase (rad)", "Phase Gradient (rad/m)")
        println("-"^110)
        res = results[k]
        for i in 1:N_z
            @printf("%6.1f | %15.6f | %20.6f | %20.6f | %20.6f\n",
                    z_tower[i], res.amplitude[i], res.phase_raw[i], res.phase_unwrapped[i], res.phase_gradient[i])
        end
        println("="^110)
        
        # Verify and classify the physical behavior of the mode
        amp = res.amplitude
        weighted_grad = sum(amp .* abs.(res.phase_gradient)) / sum(amp)
        if weighted_grad < 0.02
            @printf("Classification: >> JET STANDING MODE << (Weighted gradient = %.4f rad/m)\n", weighted_grad)
        else
            @printf("Classification: >> PROPAGATING WAVE MODE << (Weighted gradient = %.4f rad/m)\n", weighted_grad)
        end
        println("="^110)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
