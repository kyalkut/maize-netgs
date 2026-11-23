# 5. netGS step 1

Scripts (in `code/`): `netgs_step1.py`, `run_netgs_maize.py`, `merge_netgs.py`, `qp_admm.py`

- `feasibility_per_line.csv`: `max_biomass_DS` = calculated maximum growth per line (0.2842 to 0.2868),
  `dry_mass_DS_g` = measured dry mass (Zhang et al. Table S1).
- `netgs_summary_DS_paper_cvxpy_ser.json`: `s_model` (0.2855), `halfwidth` (window, 0.0030), `status` (371 solved),
  `biomass_within_window`, `median_cv_raw` (0.217), `median_cv_after_growth_normalisation` (0.0048).
- `netgs_fluxes_DS_paper_cvxpy_ser.csv`: main run, Equation 2 objective, 371 lines x 289 fluxes. `_code_` = code objective.
- `netgs_fluxes_DS_paper.csv` / `netgs_summary_DS_paper.json`: old run with our first solver (only for comparison; 67 lines invalid).
- File names: `paper` / `code` = objective, `cvxpy` = Clarabel solver, `ser` = serine pathway allowed, `co` = Rubisco ratio (section 8).
