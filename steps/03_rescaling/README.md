# 03_rescaling

log2 data -> mmol/gDW: relative factor 2^(x - WW mean), x anchor, + Table III part, x mass factor s (536.3 mg).

**Scripts:** bof_pipeline.py, steps 5-9

**Inputs:**
- `inputs/anchors.csv` (from `curation/anchors.csv`)
- `inputs/bof_tableIII.csv` (from `curation/bof_tableIII.csv`)
- `inputs/bof_components.csv` (from `output/bof_components.csv`)

**Outputs:**
- `outputs/relative_factors_WW.csv` (from `output/relative_factors_WW.csv`)
- `outputs/relative_factors_DS.csv` (from `output/relative_factors_DS.csv`)
- `outputs/free_pool_mmol_WW.csv` (from `output/free_pool_mmol_WW.csv`)
- `outputs/free_pool_mmol_DS.csv` (from `output/free_pool_mmol_DS.csv`)
- `outputs/rescaling_trace.csv` (from `output/rescaling_trace.csv`)
- `outputs/biomreaction_WW.csv` (from `output/biomreaction_WW.csv`)
- `outputs/biomreaction_DS.csv` (from `output/biomreaction_DS.csv`)
- `outputs/qc_bof_coefficients.png` (from `output/qc_bof_coefficients.png`)

Zhang log2 values: see 02_metabolite_selection/inputs.
