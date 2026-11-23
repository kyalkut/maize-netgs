# 7. Flux-sums

Scripts (in `code/`): `fluxsum.py`, `fluxsum_constraints.py`

- `fluxsum_DS_paper_cvxpy_ser.csv`: flux-sum of each measured metabolite in each line.
- `fluxsum_vs_level_DS_paper_cvxpy_ser.csv`: per metabolite `in_bof`, `n_reactions` (active reactions making it),
  `rho_level_phi` (Spearman correlation of flux-sum with measured level), `rho_level_phi_per_growth`.
- `fluxsum_constraints_summary_DS_<pfba|fba>.json`: the penalty test; `median_rho` without, with and held out;
  `n_active` (266 parsimonious, 313 plain); `status` of the problems.
