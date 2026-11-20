# netGS for maize: BOF rescaling pipeline and first netGS run

Everything needed to go from Zhang et al. (2021) metabolomics + the C4GEM maize model to
line- and condition-specific biomass reactions (BOFs), and from there to netGS step 1 (per-line fluxes).

**Start here:**
1. `MEETING_BRIEF.md`: the findings, in order, for the meeting.
2. The four figures in `output/`.
3. Everything below is reference: what each file is and how to read it.

---

## 1. Folder layout

```
bof_pipeline/
├── MEETING_BRIEF.md          findings + questions for Anika
├── README.md                 this file
├── bof_pipeline.py           MAIN PIPELINE: data -> line-specific BOFs (steps 1-10)
├── model_curation.py         fixes for the published C4GEM file (gap-fill, broken reactions, directions)
├── netgs_step1.py            Python port of netGS step 1 (1_FluxDist.m): per-line flux distributions
├── qp_admm.py                small quadratic-programming solver used by netgs_step1.py
├── run_netgs_maize.py        runs netGS step 1 on maize (DS, reference B73)
├── analysis_data.py          G x E consistency + "do BOF metabolites predict growth?"
├── figures.py                makes the meeting figures
├── curation/                 HUMAN-EDITABLE INPUTS (every decision lives here)
├── input/                    the C4GEM SBML (+ put Zhang's Additional file 2 here)
└── output/                   everything the scripts produce
```

## 2. How to run (on your computer)
1. Copy Zhang's `13059_2021_2481_MOESM2_ESM.xlsx` into `input/` and rename it `Zhang2021_AdditionalFile2.xlsx`.
2. `pip install numpy scipy openpyxl matplotlib`
3. Run, in order:
   - `python bof_pipeline.py`: BOFs. The first run reads the 70 MB xlsx (a few minutes), later runs take seconds.
   - `python analysis_data.py`, then `python figures.py`: the analyses and figures.
   - `python run_netgs_maize.py`: netGS step 1 on maize (about 30 min for 371 lines).
   - Optional: `NETGS_OBJ=code python run_netgs_maize.py` uses the objective as implemented in the original MATLAB code (see brief, finding 5).

---

## 3. The curation files (`curation/`): what you can change

| File | What it is | How to read it | When to edit |
|---|---|---|---|
| `metabolite_mapping.csv` | One row per Zhang feature that matched a model metabolite (52 rows) | `zhang_id` / `zhang_name` → `model_metabolite`; `match_type` says how it matched; `include` = yes/no | Set `include=no` to exclude a match you don't trust |
| `bof_tableIII.csv` | The C4GEM biomass equation (Dal'Molin 2010, Table III): 37 components | `coef_mmol_per_gDW` = amount per unit biomass; `mw_g_per_mol` used for the mass balance; `compartment` where the model drains it | Normally not |
| `anchors.csv` | Absolute free-pool amount per varying metabolite (31 rows) | `free_mmol_per_gDW` + `source`. Filled for 7 amino acids (Shu et al. 2022). **Empty = placeholder is used** (1 µmol/gDW) | **Main thing to complete:** add literature values with source |
| `gapfill_reactions.json` | 46 reactions added to the model | `id`, equation (`eq`), `source` (Table S1 or manual), `reversible` | Only if Anika disagrees with a direction |
| `broken_reactions.csv` | 12 reactions in the published file that create compounds from water | `equation_in_file` shows the problem (e.g. `H2O -> D-Glucose`) | No |

---

## 4. The output files (`output/`): how to read them

### Main results
**`biomreaction_WW.csv`, `biomreaction_DS.csv`: the line-specific BOFs**
- Rows = 61 BOF components; columns = `model_metabolite`, `compartment`, `varies_by_line`, then **one column per line** (385).
- Value = coefficient (mmol per unit biomass) of that component for that line.
- `varies_by_line = yes` (31 rows): built from the data, `b = bound + free_anchor × 2^(x_line − mean_WW)`.
- `varies_by_line = no`: Table III value, slightly rescaled per line so every line has the same total mass (the "fixed part scale" in the log). That's why e.g. sucrose is 0.0945 in one line and 0.0956 in another.
- *How to read:* compare a row across columns to see line differences, and the same cell between WW and DS for the drought effect. Example: B73 indole 0.00078 (WW) → 0.0155 (DS).

**`bof_components.csv`: where every BOF coefficient comes from**
- One row per component: `in_tableIII`, `varies`, `zhang_feature` (which LC-MS feature), `bound_mmol_per_gDW` (fixed part), `free_anchor_mmol_per_gDW` (anchor), `anchor_source` (literature or PLACEHOLDER), `mw`.
- *Use it to answer "where does this number come from?"*

**`netgs_fluxes_DS_paper.csv`: netGS step 1 on maize (drought)**
- Rows = reactions active in the B73 reference; `v_ref_B73` = reference flux; then one column per line (371 lines with dry mass).
- Each column is that line's flux distribution: the fluxes that satisfy its own BOF and its measured growth while staying as close as possible to B73 (paper objective).
- *These are the "intermediate traits" netGS would predict from SNPs in step 2.* Row `BIOMASS` = the line's biomass flux (follows its dry mass by construction).
- `netgs_summary_DS_paper.json`: number of lines, s_model, solver status.

### Checks
**`feasibility_per_line.csv`**: per line, the maximum biomass the curated model can reach with that line's BOF (WW, DS), next to measured DS dry mass. All > 0 means every BOF works. The values give **s_model** for the growth rescaling.

**`model_curation_log.csv`**: every change made to the published model (added / blocked / direction set) with the reason.

**`run_log.txt`**: printed log of the BOF pipeline: what was kept or dropped and why, QC results. **Read this after every run.**

**`run_summary.json`**: the same in machine-readable form.

**`sensitivity_x0.1/`, `sensitivity_x10/`**: the same pipeline with all placeholder anchors ×0.1 and ×10, to show how much the placeholders matter.

**`no24/`**: the whole pipeline (BOFs, netGS step 1 for all four setups, FluxSum, analyses) with the 24 placeholder metabolites left out (`BOF_NO_PH=1 BOF_OUT=output/no24`); BOF = Table III + 7 amino acids, 37 components. Same file names as in `output/`.

**`variant_comparison.json`, `variant_comparison_bof_DS.csv`, `variant_comparison_reactions_<tag>.csv`** (`compare_variants.py`): main run vs `no24/`. BOF part: line max biomass, spread, ranking (Spearman), shared coefficients. netGS part, per reaction: `r_lines` (correlation of the reaction's flux across lines between variants), `level_ratio_no24_over_with24`, `cv_beyond_growth_*` (CV of flux / biomass flux across lines).

### Analyses
**`analysis_GxE_bof_metabolites.csv`**: per BOF metabolite: `r_WW_DS` (does the line ranking stay the same under drought?), variance shares `pct_line` / `pct_condition` / `pct_residual_GxE`, and `drought_shift_log2`.

**`analysis_summary.json`**: prediction accuracies of DS dry mass from metabolites (BOF / random / drought-responsive / all).

### Figures
| Figure | How to read it |
|---|---|
| `qc_bof_coefficients.png` | Per varying component: spread of its coefficient across lines relative to the reference (log2). Blue = WW, orange = DS. Amino acids barely move (protein-bound part fixed); free-only metabolites carry the variation; orange shifted from blue = drought effect. |
| `fig_GxE_consistency.png` | Dots = correlation of line values WW vs DS per BOF metabolite. Left = lines reorder under drought (G×E/plasticity or noise). Dashed line = median of all 3,116 metabolites. Amino acids (orange) are lowest. |
| `fig_prediction_drymass.png` | How well different metabolite sets predict DS dry mass. The 31 BOF metabolites (orange) predict less than 31 random ones: the growth signal is mostly outside the model's scope. |
| `fig_netgs_fluxes.png` | (after the maize run) How much each reaction's flux varies across lines in netGS step 1. |

`cache_zhang.pkl`: speed-up copy of the Zhang Excel file. Ignore; delete it if you change the input file.

## 5. New in this version
- `output/rescaling_trace.csv`: every step from Zhang's log2 intensity to the final BOF coefficient (numbered columns 1–9), per condition × line × metabolite. `relative_factors_*.csv` and `free_pool_mmol_*.csv` hold columns 3 and 5 as metabolite × line tables.
- `output/netgs_fluxes_DS_paper.csv`, `netgs_summary_DS_paper.json`, `fig_netgs_fluxes.png`: netGS step 1 on all 371 drought lines (`run_netgs_maize.py`, then `merge_netgs.py`).
- `validation/`: `validate_arabidopsis.py` runs the Python port on the original netGS Arabidopsis files (set `NETGS_REPO` to the repo's `netGS/` folder, from github.com/Hao-Tong/netGS); `compare_objectives.py` compares the paper's objective with the published code's (uses the saved `arab_fluxes_*.npy`).
