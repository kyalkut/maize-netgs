# 4. From relative intensities to biomass coefficients

Scripts (in `code/`): `bof_pipeline.py (steps 5 to 9)`

- `anchors.csv`: anchor a_j per metabolite (`free_mmol_per_gDW`, `source`; empty = placeholder of 0.001 millimoles per gram).
- `molecular_weights.csv`: molecular weights M_j of the 24 new metabolites (neutral formula).
- `anchor_candidates.csv`: literature values found for 8 of the 24 (`value_nmol_per_gFW`, `source`, `location_in_source`); not used yet.
- `relative_factors_<condition>.csv`: 2^(x_ij - mean_j), lines in columns.
- `free_pool_mmol_<condition>.csv`: anchor x relative factor, before the 1 gram scaling.
- `biomreaction_<condition>.csv`: the final coefficients b_ij (rows = 61 components, columns = lines).
- `rescaling_trace.csv`: every step of the formula for every line, metabolite and condition: `3_relative_factor`,
  `4_free_anchor_mmol`, `5_free_part_mmol`, `6_bound_part_tableIII`, `7_fixed_part_scale` (= f_i), `9_final_coefficient`.
  The report's example: rows `condition` = DS, `line` = B73, `model_metabolite` = L-Tryptophan and Indole.
