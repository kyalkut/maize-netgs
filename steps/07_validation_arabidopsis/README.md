# 07_validation_arabidopsis

Python port vs published netGS on the Arabidopsis example (s_model identical, 0.0035036); paper (Eq. 2) vs code objective: same flux patterns, r = 1.00 (cvxpy).

**Scripts:** validation/validate_arabidopsis.py, validation/compare_objectives.py (NETGS_REPO = netGS repo netGS/ folder)

**Inputs:**
- none stored here

**Outputs:**
- `outputs/arab_fluxes_paper.npy` (from `validation/arab_fluxes_paper.npy`)
- `outputs/arab_fluxes_code.npy` (from `validation/arab_fluxes_code.npy`)

Inputs (biomreaction.csv, biomass_optN.csv, model) come from github.com/Hao-Tong/netGS and are not stored here.
