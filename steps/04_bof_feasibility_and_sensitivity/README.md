# 04_bof_feasibility_and_sensitivity

Max biomass of all 770 BOFs (FBA); placeholder anchors x0.1 and x10.

**Scripts:** bof_pipeline.py, step 10; sensitivity: BOF_PH_SCALE=0.1 / 10 with BOF_OUT=output/sensitivity_x0.1 / _x10

**Inputs:**
- `inputs/biomreaction_WW.csv` (from `output/biomreaction_WW.csv`)
- `inputs/biomreaction_DS.csv` (from `output/biomreaction_DS.csv`)
- `inputs/Maize_C4GEM_vs1.0.xml` (from `input/Maize_C4GEM_vs1.0.xml`)

**Outputs:**
- `outputs/feasibility_per_line.csv` (from `output/feasibility_per_line.csv`)
- `outputs/sensitivity_x0.1/` (from `output/sensitivity_x0.1`)
- `outputs/sensitivity_x10/` (from `output/sensitivity_x10`)
