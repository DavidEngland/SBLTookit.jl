### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto.jl notebook defines an interactive, reactive dashboard to analyze
# the 2D Fast-Slow Dynamical System of Turbulent Kinetic Energy (E) and Shear (S).

# ╔═╡ 11111111-1111-1111-1111-111111111111
begin
	using Pkg
	# We can let Pluto handle package installations automatically when run.
	# We use PlutoUI for sliders, OrdinaryDiffEq for stiff ODE solving, and Plots for viz.
	using PlutoUI
	using OrdinaryDiffEq
	using LinearAlgebra
	using Plots
end

# ╔═╡ 22222222-2222-2222-2222-222222222222
md"""
# 2D Fast-Slow TKE & Shear Dynamics
### Interactive Pluto.jl Dashboard for Boundary Layer Turbulence Analysis

This reactive notebook simulates the unified fast-slow system of microscale **Turbulent Kinetic Energy (TKE)** $E$ and mesoscale **Mean Wind Shear** $S$. 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000003
md"""
### 1. Interactive Model Parameters
Adjust the sliders below to dynamically bifurcate the critical manifold, vary the scale separation, and trigger shear-driven TKE collapse:
"""

# ╔═╡ 00000000-0000-0000-0000-000000000004
# Sliders
@bind G Slider(0.0:0.01:1.0, default=0.5, show_value=true, label="Geostrophic Forcing (G)")

# ╔═╡ 00000000-0000-0000-0000-000000000005
@bind epsilon Slider(0.001:0.001:0.1, default=0.01, show_value=true, label="Scale Separation (ε)")

# ╔═╡ 00000000-0000-0000-0000-000000000006
@bind cb Slider(0.1:0.05:1.0, default=0.5, show_value=true, label="Buoyant Sink Efficiency (c_b)")

# ╔═╡ 00000000-0000-0000-0000-000000000007
@bind c1 Slider(0.5:0.1:3.0, default=1.8, show_value=true, label="Shear Damping Coefficient (c_1)")

# ╔═╡ 00000000-0000-0000-0000-000000000008
@bind delta_val Slider([1e-6, 1e-5, 1e-4, 1e-3, 1e-2], default=1e-6, show_value=true, label="Regularization Floor (δ)")

# ╔═╡ 00000000-0000-0000-0000-000000000009
@bind E_init Slider(0.0:0.05:2.0, default=0.10, show_value=true, label="Initial TKE (E_0)")

# ╔═╡ 00000000-0000-0000-0000-000000000010
@bind S_init Slider(0.0:0.01:0.50, default=0.20, show_value=true, label="Initial Shear (S_0)")

# ╔═╡ 00000000-0000-0000-0000-000000000011
# Group parameters into a Dict for reactive reuse
params = begin
	Dict(
		:G => G,
		:epsilon => epsilon,
		:cb => cb,
		:c1 => c1,
		:delta => delta_val,
		:alpha => 0.1,
		:l => 0.1,
		:phi => 1.0,
		:N2 => 0.005,
		:c2 => 0.05
	)
end

# ╔═╡ 00000000-0000-0000-0000-000000000012
md"""
### 2. Dynamical Equations & Numerical Functions
The fast-slow ODE system is given by:
*   **Fast Subsystem Dynamics:**
    $$\epsilon \frac{dE}{dt} = f(E, S) = l \sqrt{E + \delta} \left(S^2 - \phi N^2 - c_b N^2 \frac{E}{E + \alpha}\right) - \frac{E^{3/2}}{l}$$
*   **Slow Subsystem Dynamics:**
    $$\frac{dS}{dt} = g(E, S) = G - c_1 E S - c_2 S$$
"""

# ╔═╡ 00000000-0000-0000-0000-000000000013
# Fast subsystem vector field f(E, S)
function f_sub(E, S, p)
	term1 = p[:l] * sqrt(E + p[:delta]) * (S^2 - p[:phi]*p[:N2] - p[:cb]*p[:N2]*E/(E + p[:alpha]))
	term2 = (E^1.5) / p[:l]
	return term1 - term2
end

# ╔═╡ 00000000-0000-0000-0000-000000000014
# Slow subsystem vector field g(E, S)
function g_sub(E, S, p)
	return p[:G] - p[:c1]*E*S - p[:c2]*S
end

# ╔═╡ 00000000-0000-0000-0000-000000000015
# Exact analytical derivative lambda_f = df/dE
function lambda_f(E, S, p)
	term1 = p[:l] * (S^2 - p[:phi]*p[:N2]) / (2.0 * sqrt(E + p[:delta]))
	term2_num = E / (2.0 * sqrt(E + p[:delta]) * (E + p[:alpha])) + (p[:alpha] * sqrt(E + p[:delta])) / ((E + p[:alpha])^2)
	term2 = p[:cb] * p[:N2] * p[:l] * term2_num
	term3 = 1.5 * sqrt(E) / p[:l]
	return term1 - term2 - term3
end

# ╔═╡ 00000000-0000-0000-0000-000000000016
# Analytical Critical Manifold S(E) where f(E, S) = 0
function get_M0_S(E, p)
	val = p[:phi]*p[:N2] + p[:cb]*p[:N2]*E/(E + p[:alpha]) + (E^1.5)/(p[:l]^2 * sqrt(E + p[:delta]))
	return val >= 0.0 ? sqrt(val) : 0.0
end

# ╔═╡ 00000000-0000-0000-0000-000000000017
# Newton-Raphson Solver to pinpoint the saddle-node fold boundary F(E, S) = [f, lambda_f]^T = 0
function find_fold(p; E0=0.1, S0=0.07, max_iter=100, tol=1e-8)
	x = [E0, S0]
	h = 1e-6
	for k in 1:max_iter
		E, S = x[1], x[2]
		E = max(E, 0.0)
		
		f_val = f_sub(E, S, p)
		lam_val = lambda_f(E, S, p)
		F = [f_val, lam_val]
		
		if norm(F) < tol
			return x, true
		end
		
		# Jacobian using finite differences
		f_E_p = f_sub(E+h, S, p)
		lam_E_p = lambda_f(E+h, S, p)
		
		f_S_p = f_sub(E, S+h, p)
		lam_S_p = lambda_f(E, S+h, p)
		
		J = [
			(f_E_p - f_val)/h   (f_S_p - f_val)/h;
			(lam_E_p - lam_val)/h (lam_S_p - lam_val)/h
		]
		
		try
			dx = J \ F
			x_next = x - dx
			x_next[1] = max(x_next[1], 0.0)
			x_next[2] = max(x_next[2], 0.0)
			
			if norm(x_next - x) < 1e-10
				x = x_next
				break
			end
			x = x_next
		catch
			break
		end
	end
	return x, false
end

# ╔═╡ 00000000-0000-0000-0000-000000000018
md"""
### 3. Stiff Numerical Integration
Using `OrdinaryDiffEq.jl`'s highly optimized stiff ODE solver `Rodas5()` to integrate the trajectory over a long time window ($t \in [0, 100]$):
"""

# ╔═╡ 00000000-0000-0000-0000-000000000019
# Set up and solve the system
sol_data = begin
	function odesystem!(du, u, p, t)
		E, S = u[1], u[2]
		du[1] = f_sub(E, S, p) / p[:epsilon]
		du[2] = g_sub(E, S, p)
	end
	
	u0 = [E_init, S_init]
	tspan = (0.0, 100.0)
	prob = ODEProblem(odesystem!, u0, tspan, params)
	solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)
end

# ╔═╡ 00000000-0000-0000-0000-000000000020
md"""
### 4. Interactive Phase Space & Stability Analysis
The plot below contains:
*   The **Critical Manifold** ($M_0$): Attracting branch (solid blue, $\lambda_f < 0$) and Repelling branch (dashed red, $\lambda_f > 0$).
*   The **Slow Nullcline** (solid green): $g(E, S) = 0$.
*   The integrated **Phase Space Trajectory** (solid black curve) with its starting point marked.
*   The exact **Saddle-Node Fold Point** $(E_{\text{fold}}, S_{\text{fold}})$ computed via the Newton-Raphson formulation (purple star).
"""

# ╔═╡ 00000000-0000-0000-0000-000000000021
# Generate Phase Portrait Plot
begin
	# Grid for Critical Manifold
	E_grid = range(0.0, 2.0, length=200)
	M0_S = [get_M0_S(E, params) for E in E_grid]
	M0_lam = [lambda_f(E, M0_S[i], params) for (i, E) in enumerate(E_grid)]
	
	# Separate into Attracting and Repelling branches
	attr_E = E_grid[M0_lam .<= 0.0]
	attr_S = M0_S[M0_lam .<= 0.0]
	
	rep_E = E_grid[M0_lam .> 0.0]
	rep_S = M0_S[M0_lam .> 0.0]
	
	# Slow Nullcline
	slow_E = range(0.001, 2.0, length=100)
	slow_S = [params[:G] / (params[:c1]*E + params[:c2]) for E in slow_E]
	
	# Trajectory
	traj_E = sol_data[1, :]
	traj_S = sol_data[2, :]
	
	# Find Fold Point
	fold_pt, fold_conv = find_fold(params)
	
	# Create Plot
	p_phase = plot(title="2D Fast-Slow Phase Space Dynamics", xlabel="Turbulent Kinetic Energy (E)", ylabel="Mean Shear (S)", xlims=(0.0, 2.0), ylims=(0.0, 0.40), legend=:topright, size=(800, 500))
	
	# Plot Attracting Branch
	plot!(p_phase, attr_E, attr_S, color=:blue, lw=2.5, label="Attracting Branch (λ_f < 0)")
	
	# Plot Repelling Branch if it exists
	if length(rep_E) > 0
		plot!(p_phase, rep_E, rep_S, color=:red, lw=2.5, linestyle=:dash, label="Repelling Branch (λ_f > 0)")
	end
	
	# Plot Slow Nullcline
	plot!(p_phase, slow_E, slow_S, color=:green, lw=2.0, label="Slow Nullcline (g = 0)")
	
	# Plot Trajectory
	plot!(p_phase, traj_E, traj_S, color=:black, lw=2.0, label="Integrated Trajectory")
	scatter!(p_phase, [E_init], [S_init], color=:orange, markersize=7, marker=:circle, label="Initial State")
	scatter!(p_phase, [traj_E[end]], [traj_S[end]], color=:cyan, markersize=7, marker=:circle, label="Steady State")
	
	# Plot Fold Point
	if fold_conv && fold_pt[1] >= 0.0
		scatter!(p_phase, [fold_pt[1]], [fold_pt[2]], color=:purple, markersize=9, marker=:star, label="Fold Point (Newton-Raphson)")
	end
	
	p_phase
end

# ╔═╡ 00000000-0000-0000-0000-000000000022
md"""
### 5. Stiff Trajectory Time Series
Below is the fast relaxation onto the manifold followed by slow drift toward the equilibrium:
"""

# ╔═╡ 00000000-0000-0000-0000-000000000023
begin
	p_time = plot(title="Time Series Evolution", xlabel="Time (t)", ylabel="Magnitude", size=(800, 350))
	plot!(p_time, sol_data.t, sol_data[1, :], color=:black, lw=2, label="TKE E(t)")
	plot!(p_time, sol_data.t, sol_data[2, :], color=:brown, lw=2, label="Shear S(t)")
	p_time
end

# ╔═╡ 00000000-0000-0000-0000-000000000024
md"""
### 6. Fast Eigenvalue Track ($\lambda_f$)
This monitors normal hyperbolicity over time. If $\lambda_f$ crosses above $0$, the manifold becomes repelling and catastrophic state transition occurs:
"""

# ╔═╡ 00000000-0000-0000-0000-000000000025
begin
	traj_lam = [lambda_f(sol_data[1, i], sol_data[2, i], params) for i in 1:length(sol_data.t)]
	p_lam = plot(title="Fast Local Eigenvalue Over Time", xlabel="Time (t)", ylabel="λ_f (df/dE)", size=(800, 300), legend=false)
	plot!(p_lam, sol_data.t, traj_lam, color=:red, lw=2)
	hline!(p_lam, [0.0], color=:black, lw=1, linestyle=:dash)
	p_lam
end

# ╔═╡ Cell order:
# ╠═ 11111111-1111-1111-1111-111111111111
# ╠═ 22222222-2222-2222-2222-222222222222
# ╠═ 00000000-0000-0000-0000-000000000003
# ╠═ 00000000-0000-0000-0000-000000000004
# ╠═ 00000000-0000-0000-0000-000000000005
# ╠═ 00000000-0000-0000-0000-000000000006
# ╠═ 00000000-0000-0000-0000-000000000007
# ╠═ 00000000-0000-0000-0000-000000000008
# ╠═ 00000000-0000-0000-0000-000000000009
# ╠═ 00000000-0000-0000-0000-000000000010
# ╠═ 00000000-0000-0000-0000-000000000011
# ╠═ 00000000-0000-0000-0000-000000000012
# ╠═ 00000000-0000-0000-0000-000000000013
# ╠═ 00000000-0000-0000-0000-000000000014
# ╠═ 00000000-0000-0000-0000-000000000015
# ╠═ 00000000-0000-0000-0000-000000000016
# ╠═ 00000000-0000-0000-0000-000000000017
# ╠═ 00000000-0000-0000-0000-000000000018
# ╠═ 00000000-0000-0000-0000-000000000019
# ╠═ 00000000-0000-0000-0000-000000000020
# ╠═ 00000000-0000-0000-0000-000000000021
# ╠═ 00000000-0000-0000-0000-000000000022
# ╠═ 00000000-0000-0000-0000-000000000023
# ╠═ 00000000-0000-0000-0000-000000000024
# ╠═ 00000000-0000-0000-0000-000000000025
