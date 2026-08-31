.PHONY: all instantiate process heatmaps gspt gabls3 paper clean test symlinks help

JULIA = julia --project=.
PIPELINE = ./run_pipeline.sh
DATA_SRC = $(abspath ../SpectralBL-Analytics/data)
LATEXMK = latexmk -pdf -interaction=nonstopmode -halt-on-error

HEATMAP_SCRIPT = scripts/plot_obukhov_heatmaps.jl
GSPT_SCRIPT = scripts/plot_gspt_transition.jl
EXPORTS_SCRIPT = scripts/run_campaign_exports.jl
GABLS3_DRIVER = scripts/gabls3-gspt-driver.jl
GABLS3_CLIM = scripts/gabls3-gspt-climatology.jl

MANUSCRIPT_DIR = paper
MANUSCRIPT_TEX = $(MANUSCRIPT_DIR)/main.tex

all: symlinks instantiate process gabls3 heatmaps gspt paper

symlinks:
	@mkdir -p data/raw
	@test -L data/datasets.json || ln -sf $(DATA_SRC)/datasets.json data/datasets.json
	@test -L data/drafts || ln -sf $(DATA_SRC)/drafts data/drafts
	@test -L data/ameriflux || ln -sf $(DATA_SRC)/ameriflux data/ameriflux
	@test -L data/raw/bllast || ln -sf $(DATA_SRC)/bllast/aeris-catalogue-prod/data/wget/113d60ba-81c9-fc79-31b2-7bbfb524fa57/processed/uniform_processing data/raw/bllast
	@test -L data/raw/floss || ln -sf $(DATA_SRC)/floss data/raw/floss
	@test -L data/raw/cases99 || ln -sf $(DATA_SRC)/cases99 data/raw/cases99
	@test -L data/raw/sheba || ln -sf $(DATA_SRC)/sheba data/raw/sheba
	@if [ -d "$(DATA_SRC)/gabls3" ]; then \
		test -L data/raw/gabls3 || ln -sf $(DATA_SRC)/gabls3 data/raw/gabls3; \
	elif [ -d "$(DATA_SRC)/gabs3" ]; then \
		test -L data/raw/gabls3 || ln -sf $(DATA_SRC)/gabs3 data/raw/gabls3; \
	fi
	@echo "Symlinks verified against $(DATA_SRC)"

instantiate:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

process:
	$(PIPELINE)

gabls3: symlinks
	@echo "Processing GABLS3 Cabauw NetCDF observations and computing GSPT coordinates..."
	$(JULIA) $(GABLS3_DRIVER)
	@echo "Compiling GABLS3 GSPT fold illusion climatology..."
	$(JULIA) $(GABLS3_CLIM)

heatmaps: symlinks
	@echo "Generating SBLToolkit L(z,t) observational heatmaps..."
	$(JULIA) $(HEATMAP_SCRIPT)

gspt: symlinks
	@echo "Generating GSPT Phase 2 dynamic (R_coord, t) transition surfaces..."
	$(JULIA) $(GSPT_SCRIPT)

exports: symlinks
	@echo "Generating campaign exports and slow manifold diagnostics..."
	$(JULIA) $(EXPORTS_SCRIPT)

paper:
	@if [ -f "$(MANUSCRIPT_TEX)" ]; then \
		echo "Compiling LaTeX manuscript..."; \
		cd $(MANUSCRIPT_DIR) && $(LATEXMK) main.tex; \
	else \
		echo "Manuscript main.tex not found in $(MANUSCRIPT_DIR)/, skipping paper build."; \
	fi

test:
	$(JULIA) -e 'using Pkg; Pkg.test()'

clean:
	rm -rf data/processed/*.jld2
	rm -rf reports/generated/sbltoolkit_heatmaps/*
	rm -rf reports/generated/gspt_phase2/*
	rm -rf reports/generated/campaign_exports/*
	rm -f gabls3_gspt_coordinates.csv gabls3_gspt_coordinates_synthetic.csv
	@if [ -d "$(MANUSCRIPT_DIR)" ]; then \
		cd $(MANUSCRIPT_DIR) && latexmk -c 2>/dev/null || true; \
	fi

help:
	@echo "SBLToolkit.jl Build Commands:"
	@echo "  make instantiate  - Install and precompile Julia dependencies"
	@echo "  make process      - Execute parallel processing pipeline across datasets"
	@echo "  make gabls3       - Ingest GABLS3 NetCDF & compile physical GSPT climatology"
	@echo "  make heatmaps     - Generate L(z,t) observational heatmaps for campaigns"
	@echo "  make gspt         - Generate dynamic GSPT (R_coord, t) transition surfaces"
	@echo "  make exports      - Run campaign data exports and manifold heatmaps"
	@echo "  make paper        - Compile LaTeX manuscript (paper/main.tex) to PDF"
	@echo "  make all          - Symlink data, process GABLS3, render figures, and build PDF"
	@echo "  make clean        - Remove processed data, generated figures, and LaTeX build files"
	@echo "  make test         - Run package test suite"
	@echo "  make symlinks     - Create necessary symlinks for raw data"