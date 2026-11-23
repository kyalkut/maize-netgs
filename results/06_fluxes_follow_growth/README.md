# 6. What step 1 gives, and why: fluxes follow growth

Scripts (in `code/`): `merge_netgs.py`, `compare_variants.py`, `bof_pipeline.py (BOF_NO_PH=1, BOF_PH_SCALE)`

- Coefficients of variation (21.7% raw, 0.48% per unit growth): `05_netgs_step1/netgs_summary_DS_paper_cvxpy_ser.json`.
  Dry mass varies by 22.2%: `feasibility_per_line.csv`, column `dry_mass_DS_g` (standard deviation / mean).
- `variant_comparison_reactions_<run>.csv`: per reaction, coefficient of variation per unit growth with and without the 24
  (76 against 8 reactions above 5%).
- `variant_comparison.json`: spread of maximum growth (0.91% against 0.17%), line order, step 1 summary of both variants.
- `variant_comparison_bof_DS.csv`: coefficients with and without the 24 (`median_with24`, `median_no24`, `spearman_lines`).
- `no24/`: the variant without the 24 placeholders (all its files are in `output/no24/` of the code folder after a run).
- `sensitivity_x0.1/`, `sensitivity_x10/`: maximum growth with the 24 placeholders 10 times smaller or larger
  (spread 0.20% and 6.89%).
