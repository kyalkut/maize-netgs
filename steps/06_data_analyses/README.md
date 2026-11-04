# 06_data_analyses

G x E consistency per BOF metabolite (WW vs DS line correlation); ridge CV prediction of DS dry mass.

**Scripts:** analysis_data.py, figures.py

**Inputs:**
- `inputs/bof_components.csv` (from `output/bof_components.csv`)
- `inputs/run_summary.json` (from `output/run_summary.json`)
- `inputs/metabolite_mapping.csv` (from `curation/metabolite_mapping.csv`)

**Outputs:**
- `outputs/analysis_GxE_bof_metabolites.csv` (from `output/analysis_GxE_bof_metabolites.csv`)
- `outputs/analysis_summary.json` (from `output/analysis_summary.json`)
- `outputs/fig_GxE_consistency.png` (from `output/fig_GxE_consistency.png`)
- `outputs/fig_prediction_drymass.png` (from `output/fig_prediction_drymass.png`)

Zhang data: see 02_metabolite_selection/inputs.
