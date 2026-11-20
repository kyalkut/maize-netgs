# 08_fluxsum

Flux-sum phi = 1/2 sum|S_ij v_j| of each measured metabolite per DS line (from step 05 fluxes); Spearman r with the measured DS log2 level across lines.

**Scripts:** fluxsum.py (NETGS_TAG picks the flux file, default paper_cvxpy_ser)

**Inputs:**
- `inputs/netgs_fluxes_DS_paper_cvxpy_ser.csv` (from `output/netgs_fluxes_DS_paper_cvxpy_ser.csv`)
- `inputs/biomreaction_DS.csv` (from `output/biomreaction_DS.csv`)
- `inputs/bof_components.csv` (from `output/bof_components.csv`)
- `inputs/metabolite_mapping.csv` (from `curation/metabolite_mapping.csv`)
- `inputs/Maize_C4GEM_vs1.0.xml` (from `input/Maize_C4GEM_vs1.0.xml`)

**Outputs:**
- `outputs/fluxsum_DS_paper_cvxpy_ser.csv` (from `output/fluxsum_DS_paper_cvxpy_ser.csv`)
- `outputs/fluxsum_vs_level_DS_paper_cvxpy_ser.csv` (from `output/fluxsum_vs_level_DS_paper_cvxpy_ser.csv`)
- `outputs/fluxsum_summary_DS_paper_cvxpy_ser.json` (from `output/fluxsum_summary_DS_paper_cvxpy_ser.json`)
- `outputs/fig_fluxsum.png` (from `output/fig_fluxsum.png`)
- `outputs/fluxsum_DS_paper_cvxpy_co.csv` (from `output/fluxsum_DS_paper_cvxpy_co.csv`)
- `outputs/fluxsum_vs_level_DS_paper_cvxpy_co.csv` (from `output/fluxsum_vs_level_DS_paper_cvxpy_co.csv`)
- `outputs/fluxsum_summary_DS_paper_cvxpy_co.json` (from `output/fluxsum_summary_DS_paper_cvxpy_co.json`)
- `outputs/fig_fluxsum_paper_cvxpy_co.png` (from `output/fig_fluxsum_paper_cvxpy_co.png`)
- `outputs/fluxsum_constraints_summary_DS_pfba.json` (from `output/fluxsum_constraints_summary_DS_pfba.json`)
- `outputs/fluxsum_constraints_fluxes_DS_pfba.csv` (from `output/fluxsum_constraints_fluxes_DS_pfba.csv`)
- `outputs/fluxsum_constraints_summary_DS_fba.json` (from `output/fluxsum_constraints_summary_DS_fba.json`)
- `outputs/fluxsum_constraints_fluxes_DS_fba.csv` (from `output/fluxsum_constraints_fluxes_DS_fba.csv`)

Zhang log2 values: see 02_metabolite_selection/inputs. Main variant excludes the biomass reaction (circular); rho_level_phi_with_bof includes it. fluxsum_constraints_*: test (fluxsum_constraints.py, NETGS_N=40 lines, FS_REF=pfba/fba) with the Table III BOF scaled to 1 g and flux-sum penalties instead of data-driven BOF components.
