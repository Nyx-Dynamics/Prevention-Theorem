# Makefile for Prevention Theorem Manuscript Figures (Science)
# Standardizes generation of Fig 1, Fig 2, Fig 3, and Supp Fig S3

PYTHON = python3
SRC_DIR = SRC
RESULTS_DIR = results/city_analysis
STOCHASTIC_RESULTS = pep_stochastic_results

.PHONY: all clean fig1 fig2 fig3 aidsvu

all: fig1 fig2 aidsvu fig3

# Script 1 — PEP_parenteral_perelson.py
# Produces: pep_parenteral_vl_sweep.png → Fig 1
fig1:
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 1..."
	cd $(SRC_DIR) && $(PYTHON) PEP_parenteral_perelson.py

# Script 2 — PEP_stochastic_perelson.py
# Produces: pep_stochastic_vl_uncertainty.png → Fig 2 + CSVs
fig2:
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 2 and Stochastic results..."
	cd $(SRC_DIR) && $(PYTHON) PEP_stochastic_perelson.py

# Script 3a — AIDSVu_city_stratified_perelson.py
# Produces: city_vl_profiles.csv, city_pep_efficacy_results.csv, etc.
aidsvu:
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Parsing AIDSVu city profiles..."
	cd $(SRC_DIR) && $(PYTHON) AIDSVu_city_stratified_perelson.py
	@cp $(RESULTS_DIR)/*.csv $(SRC_DIR)/

# Script 3b — city_stratified_figures.py
# Produces: Fig_CityStratified_PEP.png → Fig 3, Fig_CityComparison_Focus.png → Supp Fig S3
fig3: aidsvu
	@echo "[$(shell date +'%Y-%m-%d %H:%M')] Generating Fig 3 and Supp Fig S3..."
	cd $(SRC_DIR) && $(PYTHON) city_stratified_figures.py

clean:
	rm -f $(SRC_DIR)/*.png
	rm -f $(SRC_DIR)/*.csv
	rm -rf $(RESULTS_DIR)
	rm -rf $(STOCHASTIC_RESULTS)
	rm -rf SRC/pep_stochastic_results
