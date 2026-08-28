for nc_file in ../SpectralBL-Analytics/data/floss/dee002169518/datapGoGK6/*.nc; do
    echo "Processing $nc_file..."
    julia scripts/track_a_regularization.jl "$nc_file"
done