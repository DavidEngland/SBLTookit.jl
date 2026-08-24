.PHONY: all instantiate process heatmaps exports clean test symlinks help

JULIA = julia --project=.
PIPELINE = ./run_pipeline.sh
DATA_SRC = $(abspath ../SpectralBL-Analytics/data)

# Scripts
HEATMAP_SCRIPT = scripts/plot_obukhov_heatmaps.jl
EXPORTS_SCRIPT = scripts/run_campaign_exports.jl

all: symlinks instantiate process heatmaps

symlinks:
	@mkdir -p data/raw
	@test -L data/datasets.json || ln -sf $(DATA_SRC)/datasets.json data/datasets.json
	@test -L data/raw/bllast || ln -sf $(DATA_SRC)/bllast/aeris-catalogue-prod/data/wget/113d60ba-81c9-fc79-31b2-7bbfb524fa57/processed/uniform_processing data/raw/bllast
	@test -L data/raw/floss || ln -sf $(DATA_SRC)/floss data/raw/floss
	@test -L data/raw/cases99 || ln -sf $(DATA_SRC)/cases99 data/raw/cases99
	@test -L data/raw/sheba || ln -sf $(DATA_SRC)/sheba data/raw/sheba
	@test -L data/raw/gabls3 || ln -sf $(DATA_SRC)/gabls3 data/raw/gabls3
	@echo "Symlinks verified against $(DATA_SRC)"

instantiate:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

process:
	$(PIPELINE)

heatmaps:
	@echo "Generating SBLTookit $L(z,t)$ observational heatmaps..."
	$(JULIA) $(HEATMAP_SCRIPT)

exports:
	@echo "Generating campaign exports and slow manifold diagnostics..."
	$(JULIA) $(EXPORTS_SCRIPT)

test:
	$(JULIA) -e 'using Pkg; Pkg.test()'

clean:
	rm -rf data/processed/*.jld2
	rm -rf reports/generated/sbltookit_heatmaps/*
	rm -rf reports/generated/campaign_exports/*

help:
	@echo "SBLTookit.jl Build Commands:"
	@echo "  make instantiate  - Install and precompile Julia dependencies"
	@echo "  make process      - Execute parallel processing pipeline across datasets"
	@echo "  make heatmaps     - Generate L(z,t) observational heatmaps for campaigns"
	@echo "  make exports      - Run campaign data exports and manifold heatmaps"
	@echo "  make all          - Symlink data, instantiate dependencies, process, and render heatmaps"
	@echo "  make clean        - Remove processed JLD2 output files and generated figures"
	@echo "  make test         - Run package test suite"
	@echo "  make symlinks     - Create necessary symlinks for raw data"