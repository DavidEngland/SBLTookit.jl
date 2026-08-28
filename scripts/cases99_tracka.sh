for nc_file in data/raw/cases99/raw/ncar_eol_dee0099881/*.nc; do
    echo "Processing $nc_file..."
    julia scripts/track_a_regularization.jl "$nc_file"
done