for nc_file in data/raw/cases99/raw/ncar_eol_dee0099881/*.nc; do
    echo "Processing $nc_file..."
    julia src/visualize_ri_heatmaps.jl "$nc_file"
done