# 05_netgs_step1

Per DS line: s_model, growth rescaling from dry mass, B73 reference (pFBA), QP closest to reference.

**Scripts:** run_netgs_maize.py (netgs_step1.py; solver cvxpy/Clarabel, qp_admm.py = old ADMM), merge_netgs.py (figure)

**Inputs:**
- `inputs/biomreaction_DS.csv` (from `output/biomreaction_DS.csv`)
- `inputs/feasibility_per_line.csv` (from `output/feasibility_per_line.csv`)
- `inputs/Maize_C4GEM_vs1.0.xml` (from `input/Maize_C4GEM_vs1.0.xml`)

**Outputs:**
- `outputs/netgs_fluxes_DS_paper.csv` (from `output/netgs_fluxes_DS_paper.csv`)
- `outputs/netgs_summary_DS_paper.json` (from `output/netgs_summary_DS_paper.json`)
- `outputs/fig_netgs_fluxes.png` (from `output/fig_netgs_fluxes.png`)
- `output/fig_netgs_fluxes_code_cvxpy_ser.png`: not found, run the step first
- `outputs/fig_netgs_fluxes_paper_cvxpy_co.png` (from `output/fig_netgs_fluxes_paper_cvxpy_co.png`)
- `outputs/fig_netgs_fluxes_code_cvxpy_co.png` (from `output/fig_netgs_fluxes_code_cvxpy_co.png`)
- `outputs/netgs_fluxes_DS_code_cvxpy.csv` (from `output/netgs_fluxes_DS_code_cvxpy.csv`)
- `outputs/netgs_summary_DS_code_cvxpy.json` (from `output/netgs_summary_DS_code_cvxpy.json`)
- `outputs/netgs_fluxes_DS_code_cvxpy_ser.csv` (from `output/netgs_fluxes_DS_code_cvxpy_ser.csv`)
- `outputs/netgs_summary_DS_code_cvxpy_ser.json` (from `output/netgs_summary_DS_code_cvxpy_ser.json`)
- `outputs/netgs_fluxes_DS_paper_cvxpy_ser.csv` (from `output/netgs_fluxes_DS_paper_cvxpy_ser.csv`)
- `outputs/netgs_summary_DS_paper_cvxpy_ser.json` (from `output/netgs_summary_DS_paper_cvxpy_ser.json`)
- `outputs/netgs_fluxes_DS_code_cvxpy_co.csv` (from `output/netgs_fluxes_DS_code_cvxpy_co.csv`)
- `outputs/netgs_summary_DS_code_cvxpy_co.json` (from `output/netgs_summary_DS_code_cvxpy_co.json`)
- `outputs/netgs_fluxes_DS_paper_cvxpy_co.csv` (from `output/netgs_fluxes_DS_paper_cvxpy_co.csv`)
- `outputs/netgs_summary_DS_paper_cvxpy_co.json` (from `output/netgs_summary_DS_paper_cvxpy_co.json`)

Dry mass per line: column dry_mass_DS_g of feasibility_per_line.csv (from Zhang Table S1). The ADMM run (netgs_fluxes_DS_paper.csv) has 67 lines that are infeasible on the B73-active network (see the README). Files ending _ser add the phosphorylated serine pathway (R01513, R04173, R00582) to that network; all 371 lines are then feasible. code_cvxpy = objective and solver type of the original 1_FluxDist.m; paper_cvxpy = Eq. 2 of the paper (both interior point, Clarabel). Files ending _co add the maize RuBisCO bounds 22.2 <= vc/vo <= 88.0 (NETGS_CO, Wang et al. 2012) without the serine pathway; all 371 lines feasible, fluxes as _ser (per-reaction r median 1.00; 3 of 281 shared reactions below 0.9). netgs_fluxes_DS_paper.csv (ADMM) and netgs_fluxes_DS_code_cvxpy.csv (no serine pathway) are older runs with the 536.3 mg BOF, kept for comparison.
