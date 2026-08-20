#!/usr/bin/env bash
set -euo pipefail

# Configurations
RAW_DIR="./data/raw"
OUT_DIR="./data/processed"
MAX_JOBS=4  # Parallel worker processes

export JULIA_NUM_THREADS=2

mkdir -p "$OUT_DIR"

echo "=== Starting Parallel SBL Campaign Processing ==="

# Function to run julia task per file
process_file() {
    local nc_file="$1"
    local filename=$(basename "$nc_file")

    # Extract site prefix (e.g., FLOSS_2002.nc -> FLOSS)
    local site_name=$(echo "$filename" | cut -d'_' -f1 | tr '[:lower:]' '[:upper:]')

    echo "[RUNNING] Site: ${site_name} | File: ${filename}"
    julia --project=. scripts/process_campaign.jl "${site_name}" "${nc_file}" "${OUT_DIR}"
}

export -f process_file
export OUT_DIR

# Find all NetCDFs and pipe into GNU Parallel / xargs for batching
find "$RAW_DIR" -type f -name "*.nc" | xargs -I {} -P "$MAX_JOBS" bash -c 'process_file "$@"' _ {}

echo "=== Pipeline Completed Successfully ==="