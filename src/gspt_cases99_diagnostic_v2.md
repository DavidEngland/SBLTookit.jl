The v2 refactoring successfully eliminates the $\mathcal{O}(N)$ knot lookup bottleneck and guarantees safe I/O creation, though a few minor heap allocations remain inside the inner timestep loop.

**V2 Audit & Performance Breakdown**

* **Resolved Optimizations:**
* **$\mathcal{O}(\log N)$ Binary Search:** `searchsortedlast` in `evaluate_spline` eliminates linear knot search overhead during dense profile evaluations.
* **Safe Export Creation:** `mkpath(dirname(filepath))` prevents `SystemError` crashes when exporting to uncreated directories.
* **Primary Field Vector Reuse:** Moving `u_z`, `u_zz`, `theta_z`, `theta_zz`, `Ri_safe`, and `zeta` vector definitions outside the loop prevents primary field allocations across timesteps.


* **Remaining Allocation Hotspots:**
* **Inner Derivative Vectors:** `Ri_g_z = zeros(n_fine)` and `Ri_g_zz = zeros(n_fine)` (lines 202–203) are still instantiated inside the `for i_t in 1:n_t` loop, triggering garbage collection overhead on every timestep despite the comment.
* **Broadcast Array Copies:** Expressions like `C_const = R_double_prime .* (zeta_z .^ 2)` and the `CM_matrix[:, i_t]` assignment instantiate fresh temporary arrays per iteration instead of writing in-place.



**Zero-Allocation Polish**

To complete the zero-allocation loop pattern, pre-allocate `Ri_g_z`, `Ri_g_zz`, `C_const`, and `C_mapping` outside the loop and update the diagnostic assignment in-place:

```julia
    # Move remaining scratch buffers outside the 1:n_t loop
    Ri_g_z = zeros(n_fine)
    Ri_g_zz = zeros(n_fine)
    C_const = zeros(n_fine)
    C_mapping = zeros(n_fine)

    for i_t in 1:n_t
        # ... spline solver and evaluation ...

        for i_z in 1:n_fine
            Ri_g_z[i_z] = evaluate_spline(xk_ri, ak_ri, Mk_ri, hk_ri, z_fine[i_z]; order=1)
            Ri_g_zz[i_z] = evaluate_spline(xk_ri, ak_ri, Mk_ri, hk_ri, z_fine[i_z]; order=2)
        end

        # In-place curvature evaluation without heap allocation
        @. C_const = R_double_prime * (zeta_z^2)
        @. C_mapping = R_prime * zeta_zz

        K_0 = median(abs.(Ri_g_zz)) + 1e-6

        @. CM_matrix[:, i_t] = abs(C_mapping) / (abs(C_const) + abs(C_mapping) + p.epsilon_c * K_0)
    end

```