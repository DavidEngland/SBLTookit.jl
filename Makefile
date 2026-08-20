.PHONY: all instantiate process clean test symlinks help

JULIA = julia --project=.
PIPELINE = ./run_pipeline.sh
DATA_SRC = $(abspath ../SpectralBL-Analytics/data)

all: symlinks instantiate process

symlinks:
	@mkdir -p data/raw
	@test -L data/datasets.json || ln -s $(DATA_SRC)/datasets.json data/datasets.json
	@test -L data/raw/bllast || ln -s $(DATA_SRC)/bllast/aeris-catalogue-prod/data/wget/113d60ba-81c9-fc79-31b2-7bbfb524fa57/processed/uniform_processing data/raw/bllast
	@test -L data/raw/floss || ln -s $(DATA_SRC)/floss data/raw/floss
	@test -L data/raw/cases99 || ln -s $(DATA_SRC)/cases99 data/raw/cases99
	@test -L data/raw/sheba || ln -s $(DATA_SRC)/sheba data/raw/sheba
	@test -L data/raw/gabls3 || ln -s $(DATA_SRC)/gabs3 data/raw/gabls3
	@echo "Symlinks verified against $(DATA_SRC)"

instantiate:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

process:
	chmod +x $(PIPELINE)
	$(PIPELINE)

test:
	$(JULIA) -e 'using Pkg; Pkg.test()'

clean:
	rm -rf data/processed/*.jld2

help:
	@echo "SBLToolkit.jl Build Commands:"
	@echo "  make instantiate  - Install and precompile Julia dependencies"
	@echo "  make process      - Execute parallel processing pipeline across datasets"
	@echo "  make clean        - Remove processed JLD2 output files"
	@echo "  make test         - Run package test suite"
	@echo "  make symlinks     - Create necessary symlinks for raw data"