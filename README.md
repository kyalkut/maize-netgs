# netGS in maize: code and result files

This folder holds the code and the result files behind the report
"netGS in maize: what we did, where each number comes from, and why it comes out that way".
`results/` follows the report's sections; each sub-folder has a README that names the columns behind the numbers in the report.
`code/` holds every script with its curated inputs and runs as it is; `docs/` holds older explanatory notes.
`results/` is rebuilt from `code/output/` by `python make_package.py` (run in `code/`).

## Report section to files

| Report section | Result files | Scripts | Files |
|---|---|---|---|
| 1. Data and model | `results/01_data_and_model/` | `bof_pipeline.py` | 2 |
| 2. Repairing the model | `results/02_model_repair/` | `model_curation.py` | 4 |
| 3. Which measured metabolites enter the biomass reaction (351 to 31) | `results/03_metabolite_selection/` | `bof_pipeline.py` | 3 |
| 4. From relative intensities to biomass coefficients | `results/04_coefficients/` | `bof_pipeline.py` | 11 |
| 5. netGS step 1 | `results/05_netgs_step1/` | `netgs_step1.py`, `run_netgs_maize.py`, `merge_netgs.py`, `qp_admm.py` | 8 |
| 6. What step 1 gives, and why: fluxes follow growth | `results/06_fluxes_follow_growth/` | `merge_netgs.py`, `compare_variants.py`, `bof_pipeline.py` | 8 |
| 7. Flux-sums | `results/07_flux_sums/` | `fluxsum.py`, `fluxsum_constraints.py` | 8 |
| 8. A maize Rubisco ratio | `results/08_rubisco_ratio/` | `run_netgs_maize.py` | 6 |
| 9. netGS steps 2 to 4 | `results/09_netgs_steps2_to_4/` | `snp_prepare.py`, `netgs_steps234.py`, `analysis_data.py` | 46 |
| 10. Other maize models (open questions) | `results/10_other_models/` | `compare_models.py` | 2 |
| Validation on Arabidopsis (sections 5 and 6) | `results/11_validation_arabidopsis/` | `validation/validate_arabidopsis.py`, `validation/compare_objectives.py` | 0 |

## Scripts

| Script | What it does |
|---|---|
| `bof_pipeline.py` | Data to line-specific biomass reactions (steps 1 to 10); thresholds in PARAM at the top. Calls model_curation.py. |
| `model_curation.py` | Repairs of C4GEM: added, switched-off and one-directional reactions. |
| `netgs_step1.py` | Python port of netGS 1_FluxDist.m (maximum growth, s_model, growth target, reference, quadratic program). |
| `run_netgs_maize.py` | Runs netGS step 1 on the drought lines (settings by environment variables, see README). |
| `merge_netgs.py` | Merges step 1 runs, computes the coefficients of variation (raw and per unit growth), figure. |
| `qp_admm.py` | Our first solver (alternating direction method of multipliers); replaced by Clarabel, kept for comparison. |
| `fluxsum.py` | Flux-sums of the measured metabolites and their correlation with the measured levels. |
| `fluxsum_constraints.py` | Test: flux-sum penalty instead of data-driven biomass components (40 lines). |
| `compare_variants.py` | Main biomass reaction against the variant without the 24 placeholders. |
| `snp_prepare.py` | Markers (ZEAMAP): filtering, line matching, relationship matrix K. |
| `netgs_steps234.py` | Python port of netGS 2_rrBLUP.R, 3_Biomass.m and 4_Correlation.R; accuracy against classical genomic selection. |
| `analysis_data.py` | Metabolites predicting dry mass (ridge regression); well-watered against drought consistency. |
| `compare_models.py` | C4GEM against Simons et al. 2014 and Bogart and Myers 2016. |
| `figures.py` | Figures. |
| `validation/validate_arabidopsis.py` | Reproduces netGS step 1 on the authors' Arabidopsis data (needs the netGS repository). |
| `validation/compare_objectives.py` | Equation 2 against the code objective on Arabidopsis. |
| `validation/sbml.py` | Small reader for SBML model files. |

## Editable inputs (`code/curation/`)

| File | Content |
|---|---|
| `metabolite_mapping.csv` | Zhang metabolite name to C4GEM name |
| `anchors.csv` | anchor per metabolite (millimoles per gram dry weight) and its source |
| `bof_tableIII.csv` | C4GEM biomass reaction from Dal'Molin et al. 2010, Table III |
| `molecular_weights.csv` | molecular weights of the 24 new metabolites |
| `gapfill_reactions.json` | reactions added to the model |
| `broken_reactions.csv` | reactions switched off |

## How to run (from `code/`)

```
pip install numpy scipy openpyxl matplotlib cvxpy clarabel
python bof_pipeline.py            # sections 2-4: model repair, selection, biomass reactions (first run reads the xlsx, a few minutes)
python run_netgs_maize.py         # section 5: netGS step 1 (about 30 minutes; NETGS_OBJ=paper|code, NETGS_SOLVER=cvxpy)
python merge_netgs.py             # section 6: coefficients of variation, figure
python fluxsum.py                 # section 7
NETGS_CO=22.2:88.0 NETGS_EXTRA=none python run_netgs_maize.py   # section 8: Rubisco ratio
BOF_NO_PH=1 BOF_OUT=output/no24 python bof_pipeline.py          # section 6: variant without the 24 (then step 1 with the same BOF_OUT)
python compare_variants.py        # section 6
python snp_prepare.py             # section 9: needs the ZEAMAP marker file in input/snp/
python netgs_steps234.py          # section 9
python analysis_data.py; python compare_models.py; python figures.py
```

Results are written to `code/output/`. Inputs not included (size or third-party licence): Zhang et al. Additional file 2
(`code/input/Zhang2021_AdditionalFile2.xlsx`, doi:10.1186/s13059-021-02481-1), the ZEAMAP marker file (`code/input/snp/`),
the other maize models (`code/input/other_models/`).
