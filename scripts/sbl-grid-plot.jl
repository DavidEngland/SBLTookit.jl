# ==============================================================================
# SBL GSPT Grid-Stretching Plotting Script
# Uses Plots.jl (GR backend) to visualize the time-evolution of the adaptive 
# vertical physical grid over a descending Low-Level Jet (LLJ) core.
# ==============================================================================

using Plots
using DelimitedFiles

# Set up publication-quality plotting themes and configurations
theme(:default)
default(
    fontfamily = "sans-serif",
    tickfontsize = 10,
    guidefontsize = 11,
    titlefontsize = 12,
    legendfontsize = 9,
    linewidth = 1.5,
    grid = true,
    gridalpha = 0.3,
    gridstyle = :dash
)

function plot_adaptive_grid_evolution(csv_path::String; output_path::String="sbl-grid-stretching.png")
    # 1. Read the exported grid trajectory matrix
    # Format is N_levels x timesteps (e.g., 38 rows by 8 columns)
    if !isfile(csv_path)
        error("Grid history CSV file not found at path: $csv_path")
    end
    grid_data = readdlm(csv_path, ',', Float64)
    N_levels, timesteps = size(grid_data)
    
    # 2. Define coordinate domains
    t_steps = collect(1:timesteps)
    z_min, z_max = grid_data[1, 1], grid_data[end, 1]
    
    # 3. Simulate/compute the descending jet nose height over time
    # This matches the physical forcing model in grid-stretching-engine.jl
    llj_nose = [140.0 - (t - 1) * 14.0 for t in t_steps]
    
    # 4. Initialize Plot
    # We use a custom size and clean DPI for publication quality
    p = plot(
        size=(850, 600),
        dpi=300,
        title="Dynamic SCM Grid Stretching Over Descending Low-Level Jet Nose",
        xlabel="Time of Nocturnal Cycle (hours)",
        ylabel="Height Above Ground z (m)",
        xlims=(0.8, timesteps + 0.2),
        ylims=(z_min, z_max),
        yscale=:log10, # Logarithmic vertical scale aligns with near-surface SBL physics
        legend=:topright,
        bg=:white,
        fg=:black
    )
    
    # 5. Plot the time-evolution of the physical grid lines
    # To show the grid contraction clearly, we draw each level as a continuous time-series line
    for i in 1:N_levels
        # Use subtle gray/blue horizontal grid lines to represent vertical levels over time
        label_str = (i == 1 || i == N_levels) ? (i == 1 ? "Boundary Nodes" : "") : ""
        plot!(p, t_steps, grid_data[i, :], 
              color=:dodgerblue, 
              alpha=0.4, 
              linewidth=1.0, 
              marker=:circle, 
              markersize=2, 
              markeralpha=0.6,
              label=label_str)
    end
    
    # 6. Overlay the physical forcing trajectory (Descending LLJ Nose)
    plot!(p, t_steps, llj_nose, 
          color=:black, 
          linewidth=3.0, 
          line=:solid, 
          label="Low-Level Jet (LLJ) Nose (U_z = 0)",
          annotation=(timesteps / 2, 50.0, text("Grid lines contract dynamically\naround the vanishing-shear core", 9, :darkgray, :center)))
          
    # 7. Add visual markers indicating high-density coordinate compression
    for (t, z_n) in zip(t_steps, llj_nose)
        scatter!(p, [t], [z_n], 
                 color=:crimson, 
                 markersize=5, 
                 label= t == 1 ? "Highest Compression Zone (ζ_zz)" : "")
    end
    
    # 8. Customize Y-axis formatting for logarithmic scale readability
    yticks_vals = [2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 100.0, 140.0, 200.0]
    yticks!(p, yticks_vals, string.(yticks_vals))
    
    # Save the generated visualization to the specified outbox path
    savefig(p, output_path)
    println("[SBL PLOTTER] Successfully generated and saved: $output_path")
    return p
end

# To execute the script in a standard Julia terminal:
# plot_adaptive_grid_evolution("grid_history.csv")
