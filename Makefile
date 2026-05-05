# Makefile for Prevention Theorem (Science Advances v2 Revisions)
# Author: Junie (Autonomous Programmer)
# Date: 2026-05-05

PYTHON = python3
REPRO_SCRIPT = reproduce_all_v2.py
FIGURES_DIR = v2_revision/figures
RUNS_DIR = v2_revision/runs
SRC_DIR = SRC

.PHONY: all figure1 figure2 figure3 clean clean-runs help

all: ## Reproduce all figures (1, 2, 3, 4) into a new timestamped run directory
	$(PYTHON) $(REPRO_SCRIPT)

figure1: ## Generate Figure 1 only (Route comparison)
	$(PYTHON) v2_revision/make_figure_1_v2.py

figure2: ## Generate Figure 2 only (Stochastic VL)
	PYTHONPATH=$(SRC_DIR) $(PYTHON) $(SRC_DIR)/perelson/stochastic/PEP_stochastic_perelson.py
	@if [ -f pep_stochastic_vl_uncertainty.png ]; then \
		mv pep_stochastic_vl_uncertainty.png $(FIGURES_DIR)/Figure_2_Stochastic_Efficacy_VL_Uncertainty.png; \
		echo "Moved to $(FIGURES_DIR)/Figure_2_Stochastic_Efficacy_VL_Uncertainty.png"; \
	fi

figure3: ## Generate Figure 3 & 4 only (City stratified)
	$(PYTHON) v2_revision/city_stratified_figures.py
	@if [ -f results/city_analysis/Fig_CityStratified_PEP.png ]; then \
		mv results/city_analysis/Fig_CityStratified_PEP.png $(FIGURES_DIR)/Figure_3_City_Stratified_PEP_Efficacy.png; \
	fi
	@if [ -f results/city_analysis/Fig_CityComparison_Focus.png ]; then \
		mv results/city_analysis/Fig_CityComparison_Focus.png $(FIGURES_DIR)/Figure_S1_CityComparison_Focus.png; \
	fi

clean: ## Remove generated figures from static figures/ directory
	rm -rf $(FIGURES_DIR)/*.png
	rm -rf $(FIGURES_DIR)/*.pdf
	rm -rf pep_stochastic_results/*.csv
	rm -f pep_stochastic_vl_uncertainty.png
	@echo "Cleanup of static figures complete."

clean-runs: ## Remove all timestamped run directories
	rm -rf $(RUNS_DIR)/*
	@echo "Cleanup of all runs complete."

help: ## Display this help message
	@grep -E '^[a-zA-Z1-9_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'
