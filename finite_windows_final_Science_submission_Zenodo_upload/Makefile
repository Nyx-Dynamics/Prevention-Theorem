# Makefile for Prevention Theorem Manuscript Figures (Science)
# Standardizes generation of Fig 1, Fig 2, Fig 3, and Supp Fig S3
# Configured for flat Zenodo repository structure

PYTHON = python3
VENV = .venv
BIN = $(VENV)/bin
VENV_PYTHON = $(BIN)/python3
VENV_PIP = $(BIN)/pip

# Use VENV_PYTHON if it exists, otherwise fallback to system PYTHON
PYTHON_CMD = $(if $(wildcard $(VENV_PYTHON)),$(VENV_PYTHON),$(PYTHON))

SRC_DIR = .
RESULTS_DIR = results/city_analysis
STOCHASTIC_RESULTS = pep_stochastic_results
TIMESTAMP := $(shell date +'%Y-%m-%d_%H%M')
OUTPUT_BASE = outputs
RUN_DIR = $(OUTPUT_BASE)/$(TIMESTAMP)

.PHONY: all clean fig1 fig2 fig3 aidsvu init_run venv clean-venv

all: venv init_run fig1 fig2 aidsvu fig3
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Run complete. All outputs saved to $(RUN_DIR)"

init_run: venv
	@mkdir -p $(RUN_DIR)/$(RESULTS_DIR)
	@mkdir -p $(RUN_DIR)/$(STOCHASTIC_RESULTS)

venv: $(VENV)/bin/activate

$(VENV)/bin/activate: requirements.txt
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Creating virtual environment..."
	$(PYTHON) -m venv $(VENV)
	$(VENV_PIP) install --upgrade pip
	$(VENV_PIP) install -r requirements.txt
	@touch $(VENV)/bin/activate

# Script 1 — PEP_parenteral_perelson.py
# Produces: pep_parenteral_vl_sweep.png → Fig 1
fig1: init_run
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 1..."
	$(PYTHON_CMD) PEP_parenteral_perelson.py
	@mv FIgure_1_pep_parenteral_vl_sweep.png $(RUN_DIR)/ 2>/dev/null || true

# Script 2 — PEP_stochastic_perelson.py
# Produces: pep_stochastic_vl_uncertainty.png → Fig 2 + CSVs
fig2: init_run
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 2 and Stochastic results..."
	$(PYTHON_CMD) PEP_stochastic_perelson.py
	@mv Figure_2_pep_stochastic_vl_uncertainty.png $(RUN_DIR)/ 2>/dev/null || true
	@if [ -d "$(STOCHASTIC_RESULTS)" ]; then cp -r $(STOCHASTIC_RESULTS)/* $(RUN_DIR)/$(STOCHASTIC_RESULTS)/; fi

# Script 3a — aidsvu_city_profiles.py
# Produces: city_vl_profiles.csv, city_pep_efficacy_results.csv, etc.
aidsvu: init_run
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Parsing AIDSVu city profiles..."
	$(PYTHON_CMD) aidsvu_city_profiles.py
	@if [ -d "$(RESULTS_DIR)" ]; then cp $(RESULTS_DIR)/*.csv $(RUN_DIR)/$(RESULTS_DIR)/; fi
	@cp *.csv $(RUN_DIR)/ 2>/dev/null || true

# Script 3b — city_stratified_figures.py
# Produces: Fig_CityStratified_PEP.png → Fig 3, Fig_CityComparison_Focus.png → Supp Fig S3
fig3: aidsvu
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 3 and Supp Fig S3..."
	$(PYTHON_CMD) city_stratified_figures.py
	@if [ -d "$(RESULTS_DIR)" ]; then cp $(RESULTS_DIR)/* $(RUN_DIR)/$(RESULTS_DIR)/; fi

clean:
	rm -f *.png
	rm -f *.csv
	rm -rf $(RESULTS_DIR)
	rm -rf $(STOCHASTIC_RESULTS)
	rm -rf pep_stochastic_results
	@echo "Cleaned workspace. Run history in $(OUTPUT_BASE)/ is preserved."

clean-venv:
	rm -rf $(VENV)
	@echo "Removed virtual environment."
