.PHONY: all instantiate process clean test help

JULIA = julia --project=.
PIPELINE = ./run_pipeline.sh

all: instantiate process

instantiate:
	$(JULIA) -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

process:
	chmod +x $(PIPELINE)$(PIPELINE)

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