# 3. Which measured metabolites enter the biomass reaction (351 to 31)

Scripts (in `code/`): `bof_pipeline.py (steps 2 to 4)`

- `metabolite_mapping.csv`: Zhang name to C4GEM name (52 features, 50 metabolites).
- `run_log.txt`, steps 3 and 4: producible metabolites (35) and the data filters (31 kept); which ones were dropped and why.
- `bof_components.csv`: the 61 components of every biomass reaction. Column `in_tableIII` (yes for the 37 Table III components),
  `varies` (yes for the 31 line-specific ones), `zhang_feature`, `bound_mmol_per_gDW` (Table III part),
  `free_anchor_mmol_per_gDW` (anchor), `anchor_source`, `mw` (molecular weight).
