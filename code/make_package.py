"""Rebuild ../results/ and ../README.md from this folder (code/), following the sections of the combined report
"netGS in maize: what we did, where each number comes from, and why it comes out that way".

  README.md            overview, how to run, map: report section -> scripts -> files
  code/                every script + curation inputs + C4GEM model (writes code/output/, not tracked by git)
  results/NN_section/  copies of the result files each report section uses, with a README naming the columns

Rerun after any pipeline run:  python make_package.py  (from code/)"""
import os, shutil, glob

HERE = os.path.dirname(os.path.abspath(__file__))
DST = os.path.dirname(HERE)   # repository root

# ---------------------------------------------------------------- code
SCRIPTS = [
    ('bof_pipeline.py', 'Data to line-specific biomass reactions (steps 1 to 10); thresholds in PARAM at the top. Calls model_curation.py.'),
    ('model_curation.py', 'Repairs of C4GEM: added, switched-off and one-directional reactions.'),
    ('netgs_step1.py', 'Python port of netGS 1_FluxDist.m (maximum growth, s_model, growth target, reference, quadratic program).'),
    ('run_netgs_maize.py', 'Runs netGS step 1 on the drought lines (settings by environment variables, see README).'),
    ('merge_netgs.py', 'Merges step 1 runs, computes the coefficients of variation (raw and per unit growth), figure.'),
    ('qp_admm.py', 'Our first solver (alternating direction method of multipliers); replaced by Clarabel, kept for comparison.'),
    ('fluxsum.py', 'Flux-sums of the measured metabolites and their correlation with the measured levels.'),
    ('fluxsum_constraints.py', 'Test: flux-sum penalty instead of data-driven biomass components (40 lines).'),
    ('compare_variants.py', 'Main biomass reaction against the variant without the 24 placeholders.'),
    ('snp_prepare.py', 'Markers (ZEAMAP): filtering, line matching, relationship matrix K.'),
    ('netgs_steps234.py', 'Python port of netGS 2_rrBLUP.R, 3_Biomass.m and 4_Correlation.R; accuracy against classical genomic selection.'),
    ('analysis_data.py', 'Metabolites predicting dry mass (ridge regression); well-watered against drought consistency.'),
    ('compare_models.py', 'C4GEM against Simons et al. 2014 and Bogart and Myers 2016.'),
    ('figures.py', 'Figures.'),
    ('validation/validate_arabidopsis.py', 'Reproduces netGS step 1 on the authors\' Arabidopsis data (needs the netGS repository).'),
    ('validation/compare_objectives.py', 'Equation 2 against the code objective on Arabidopsis.'),
    ('validation/sbml.py', 'Small reader for SBML model files.'),
]
CURATION = ['metabolite_mapping.csv', 'anchors.csv', 'bof_tableIII.csv', 'molecular_weights.csv',
            'gapfill_reactions.json', 'broken_reactions.csv']

# ---------------------------------------------------------------- results per report section
O = 'output/'
SECTIONS = [
    dict(name='01_data_and_model', title='1. Data and model',
         scripts=['bof_pipeline.py'],
         files=['input/Maize_C4GEM_vs1.0.xml', 'curation/bof_tableIII.csv'],
         readme="""C4GEM model file (Dal'Molin et al. 2010, Supplemental File S1) and its biomass reaction, which exists only in the
paper's Table III (`bof_tableIII.csv`: 37 components, millimoles per gram dry weight; maintenance 30 millimoles adenosine triphosphate).
Zhang et al. (2021) Additional file 2 (72 MB) is not copied; download it from Genome Biology (doi:10.1186/s13059-021-02481-1)
into `code/input/Zhang2021_AdditionalFile2.xlsx`. Dry mass per drought line comes from its Table S1."""),
    dict(name='02_model_repair', title='2. Repairing the model',
         scripts=['model_curation.py'],
         files=['curation/gapfill_reactions.json', 'curation/broken_reactions.csv', O + 'model_curation_log.csv', O + 'run_log.txt'],
         readme="""- `gapfill_reactions.json`: the 46 added reactions (43 from Dal'Molin et al. Table S1, plus 3).
- `broken_reactions.csv`: the 12 reactions switched off, with the reason for each.
- `model_curation_log.csv` (columns `reaction`, `change`): every change made, including the 28 directions set by enzyme class.
- `run_log.txt`, step 1: counts and the checks (no adenosine triphosphate or biomass without light or carbon dioxide)."""),
    dict(name='03_metabolite_selection', title='3. Which measured metabolites enter the biomass reaction (351 to 31)',
         scripts=['bof_pipeline.py (steps 2 to 4)'],
         files=['curation/metabolite_mapping.csv', O + 'bof_components.csv', O + 'run_log.txt'],
         readme="""- `metabolite_mapping.csv`: Zhang name to C4GEM name (52 features, 50 metabolites).
- `run_log.txt`, steps 3 and 4: producible metabolites (35) and the data filters (31 kept); which ones were dropped and why.
- `bof_components.csv`: the 61 components of every biomass reaction. Column `in_tableIII` (yes for the 37 Table III components),
  `varies` (yes for the 31 line-specific ones), `zhang_feature`, `bound_mmol_per_gDW` (Table III part),
  `free_anchor_mmol_per_gDW` (anchor), `anchor_source`, `mw` (molecular weight)."""),
    dict(name='04_coefficients', title='4. From relative intensities to biomass coefficients',
         scripts=['bof_pipeline.py (steps 5 to 9)'],
         files=['curation/anchors.csv', 'curation/molecular_weights.csv', 'literature/anchor_sources/anchor_candidates.csv',
                O + 'relative_factors_WW.csv', O + 'relative_factors_DS.csv', O + 'free_pool_mmol_WW.csv', O + 'free_pool_mmol_DS.csv',
                O + 'biomreaction_WW.csv', O + 'biomreaction_DS.csv', O + 'rescaling_trace.csv', O + 'qc_bof_coefficients.png'],
         readme="""- `anchors.csv`: anchor a_j per metabolite (`free_mmol_per_gDW`, `source`; empty = placeholder of 0.001 millimoles per gram).
- `molecular_weights.csv`: molecular weights M_j of the 24 new metabolites (neutral formula).
- `anchor_candidates.csv`: literature values found for 8 of the 24 (`value_nmol_per_gFW`, `source`, `location_in_source`); not used yet.
- `relative_factors_<condition>.csv`: 2^(x_ij - mean_j), lines in columns.
- `free_pool_mmol_<condition>.csv`: anchor x relative factor, before the 1 gram scaling.
- `biomreaction_<condition>.csv`: the final coefficients b_ij (rows = 61 components, columns = lines).
- `rescaling_trace.csv`: every step of the formula for every line, metabolite and condition: `3_relative_factor`,
  `4_free_anchor_mmol`, `5_free_part_mmol`, `6_bound_part_tableIII`, `7_fixed_part_scale` (= f_i), `9_final_coefficient`.
  The report's example: rows `condition` = DS, `line` = B73, `model_metabolite` = L-Tryptophan and Indole."""),
    dict(name='05_netgs_step1', title='5. netGS step 1',
         scripts=['netgs_step1.py', 'run_netgs_maize.py', 'merge_netgs.py', 'qp_admm.py'],
         files=[O + 'feasibility_per_line.csv', O + 'netgs_fluxes_DS_paper_cvxpy_ser.csv', O + 'netgs_fluxes_DS_code_cvxpy_ser.csv',
                O + 'netgs_summary_DS_paper_cvxpy_ser.json', O + 'netgs_summary_DS_code_cvxpy_ser.json', O + 'fig_netgs_fluxes.png',
                O + 'netgs_fluxes_DS_paper.csv', O + 'netgs_summary_DS_paper.json'],
         readme="""- `feasibility_per_line.csv`: `max_biomass_DS` = calculated maximum growth per line (0.2842 to 0.2868),
  `dry_mass_DS_g` = measured dry mass (Zhang et al. Table S1).
- `netgs_summary_DS_paper_cvxpy_ser.json`: `s_model` (0.2855), `halfwidth` (window, 0.0030), `status` (371 solved),
  `biomass_within_window`, `median_cv_raw` (0.217), `median_cv_after_growth_normalisation` (0.0048).
- `netgs_fluxes_DS_paper_cvxpy_ser.csv`: main run, Equation 2 objective, 371 lines x 289 fluxes. `_code_` = code objective.
- `netgs_fluxes_DS_paper.csv` / `netgs_summary_DS_paper.json`: old run with our first solver (only for comparison; 67 lines invalid).
- File names: `paper` / `code` = objective, `cvxpy` = Clarabel solver, `ser` = serine pathway allowed, `co` = Rubisco ratio (section 8)."""),
    dict(name='06_fluxes_follow_growth', title='6. What step 1 gives, and why: fluxes follow growth',
         scripts=['merge_netgs.py', 'compare_variants.py', 'bof_pipeline.py (BOF_NO_PH=1, BOF_PH_SCALE)'],
         files=[O + 'variant_comparison.json', O + 'variant_comparison_bof_DS.csv',
                O + 'variant_comparison_reactions_paper_cvxpy_ser.csv', O + 'variant_comparison_reactions_code_cvxpy_ser.csv',
                O + 'no24/feasibility_per_line.csv', O + 'no24/netgs_summary_DS_paper_cvxpy_ser.json',
                O + 'sensitivity_x0.1/feasibility_per_line.csv', O + 'sensitivity_x10/feasibility_per_line.csv'],
         readme="""- Coefficients of variation (21.7% raw, 0.48% per unit growth): `05_netgs_step1/netgs_summary_DS_paper_cvxpy_ser.json`.
  Dry mass varies by 22.2%: `feasibility_per_line.csv`, column `dry_mass_DS_g` (standard deviation / mean).
- `variant_comparison_reactions_<run>.csv`: per reaction, coefficient of variation per unit growth with and without the 24
  (76 against 8 reactions above 5%).
- `variant_comparison.json`: spread of maximum growth (0.91% against 0.17%), line order, step 1 summary of both variants.
- `variant_comparison_bof_DS.csv`: coefficients with and without the 24 (`median_with24`, `median_no24`, `spearman_lines`).
- `no24/`: the variant without the 24 placeholders (all its files are in `output/no24/` of the code folder after a run).
- `sensitivity_x0.1/`, `sensitivity_x10/`: maximum growth with the 24 placeholders 10 times smaller or larger
  (spread 0.20% and 6.89%)."""),
    dict(name='07_flux_sums', title='7. Flux-sums',
         scripts=['fluxsum.py', 'fluxsum_constraints.py'],
         files=[O + 'fluxsum_DS_paper_cvxpy_ser.csv', O + 'fluxsum_vs_level_DS_paper_cvxpy_ser.csv', O + 'fluxsum_summary_DS_paper_cvxpy_ser.json',
                O + 'fig_fluxsum.png', O + 'fluxsum_constraints_summary_DS_pfba.json', O + 'fluxsum_constraints_summary_DS_fba.json',
                O + 'fluxsum_constraints_fluxes_DS_pfba.csv', O + 'fluxsum_constraints_fluxes_DS_fba.csv'],
         readme="""- `fluxsum_DS_paper_cvxpy_ser.csv`: flux-sum of each measured metabolite in each line.
- `fluxsum_vs_level_DS_paper_cvxpy_ser.csv`: per metabolite `in_bof`, `n_reactions` (active reactions making it),
  `rho_level_phi` (Spearman correlation of flux-sum with measured level), `rho_level_phi_per_growth`.
- `fluxsum_constraints_summary_DS_<pfba|fba>.json`: the penalty test; `median_rho` without, with and held out;
  `n_active` (266 parsimonious, 313 plain); `status` of the problems."""),
    dict(name='08_rubisco_ratio', title='8. A maize Rubisco ratio',
         scripts=['run_netgs_maize.py (NETGS_CO=22.2:88.0, NETGS_EXTRA=none)'],
         files=[O + 'netgs_fluxes_DS_paper_cvxpy_co.csv', O + 'netgs_fluxes_DS_code_cvxpy_co.csv',
                O + 'netgs_summary_DS_paper_cvxpy_co.json', O + 'netgs_summary_DS_code_cvxpy_co.json',
                O + 'fig_netgs_fluxes_paper_cvxpy_co.png', O + 'fluxsum_vs_level_DS_paper_cvxpy_co.csv'],
         readme="""- `netgs_summary_DS_paper_cvxpy_co.json`: `s_model` 0.2810, `status` (371 solved without the serine pathway).
- `netgs_fluxes_DS_<objective>_cvxpy_co.csv`: fluxes with the Rubisco ratio; compare with `05_netgs_step1` (same columns).
- Ratio values: Wang et al. 2012, Table 2 (22.2 to 88.0)."""),
    dict(name='09_netgs_steps2_to_4', title='9. netGS steps 2 to 4',
         scripts=['snp_prepare.py', 'netgs_steps234.py', 'analysis_data.py'],
         files=[O + 'snp_summary.json', O + 'snp_line_match.csv', O + 'netgs_foldid_DS.csv',
                O + 'netgs_step2_accuracy_DS_paper_cvxpy_ser.csv']
               + [O + os.path.basename(p) for p in sorted(glob.glob(os.path.join(HERE, 'output', 'netgs_step34_*.json')))]
               + [O + os.path.basename(p) for p in sorted(glob.glob(os.path.join(HERE, 'output', 'netgs_step34_folds_*.csv')))]
               + [O + 'no24/' + os.path.basename(p) for p in sorted(glob.glob(os.path.join(HERE, 'output', 'no24', 'netgs_step34_*.json')))]
               + [O + 'sensitivity_x0.1/netgs_step34_DS_paper_cvxpy_ser_paper.json', O + 'sensitivity_x10/netgs_step34_DS_paper_cvxpy_ser_paper.json',
                  O + 'analysis_summary.json', O + 'fig_prediction_drymass.png'],
         readme="""- `snp_summary.json`, `snp_line_match.csv`: markers kept (796,846) and how the 353 lines were matched (`match`).
  The marker file itself (208 MB) is not copied: ZEAMAP, ftp.cngb.org/pub/CNSA/data3/CNP0001565/zeamap/02_Variants/SNPs/AMP_SNP_anno.vcf.gz.
- `netgs_foldid_DS.csv`: the 50 x 3-fold split, the same for every run.
- `netgs_step2_accuracy_DS_paper_cvxpy_ser.csv`: per flux `accuracy_mean` (step 2), `r_with_biomass_flux`.
- `netgs_step34_DS_<step 1 run>_<step 3 objective>.json`: `netgs_accuracy_mean`, `classical_gs_drymass_mean` (0.325),
  `share_folds_netgs_better`. Suffix `_sel-<mode>`: which reactions enter step 3 (biomass, growth, all, inner, nongrowth).
  `no24/`: without the 24 placeholders. `sensitivity_x*/`: placeholders 10 times smaller or larger.
- `netgs_step34_folds_*.csv`: the same per test set (`r_netgs`, `r_classical_gs_drymass`).
- `analysis_summary.json`, `prediction`: metabolites predicting dry mass (all 0.54, 31 random 0.29, our 31 0.19)."""),
    dict(name='10_other_models', title='10. Other maize models (open questions)',
         scripts=['compare_models.py'],
         files=[O + 'model_comparison_summary.json', O + 'model_comparison_zhang_overlap.csv'],
         readme="""- `model_comparison_summary.json`: reactions, cell types, genes and biomass reaction of C4GEM, Simons et al. 2014,
  iEB5204 and iEB2140x2 (Bogart and Myers 2016).
- `model_comparison_zhang_overlap.csv`: which of Zhang's named metabolites each model contains (same automatic name rule).
  The other models' files are not copied (publishers' supplements; see `code/input/other_models/README.md`)."""),
    dict(name='11_validation_arabidopsis', title='Validation on Arabidopsis (sections 5 and 6)',
         scripts=['validation/validate_arabidopsis.py', 'validation/compare_objectives.py'],
         files=[],
         readme="""Reproduces s_model = 0.0035036 of Tong et al. (2020) on their Arabidopsis data and compares Equation 2 with
the code objective (median per-reaction correlation 1.00). Needs the authors' repository:
`git clone https://github.com/Hao-Tong/netGS`, then `NETGS_REPO=netGS/netGS/ python validation/validate_arabidopsis.py` from `code/`."""),
]

RUN = """## How to run (from `code/`)

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
"""


def cp(src_rel, dst_dir):
    src = os.path.join(HERE, src_rel)
    if not os.path.exists(src):
        print('  missing:', src_rel); return False
    os.makedirs(dst_dir, exist_ok=True)
    d = os.path.join(dst_dir, os.path.basename(src_rel.rstrip('/')))
    # keep sub-folder name for variant files (no24/, sensitivity_x*/)
    parts = src_rel.split('/')
    if len(parts) >= 3 and parts[0] == 'output':
        d = os.path.join(dst_dir, parts[1], parts[-1]); os.makedirs(os.path.dirname(d), exist_ok=True)
    shutil.copy2(src, d); return True


def main():
    res = os.path.join(DST, 'results')
    if os.path.exists(res): shutil.rmtree(res)

    rows = []
    for sec in SECTIONS:
        d = os.path.join(DST, 'results', sec['name']); os.makedirs(d, exist_ok=True)
        n = sum(cp(f, d) for f in sec['files'])
        with open(os.path.join(d, 'README.md'), 'w', encoding='utf-8') as fh:
            fh.write(f"# {sec['title']}\n\nScripts (in `code/`): {', '.join('`' + s + '`' for s in sec['scripts'])}\n\n{sec['readme']}\n")
        rows.append(f"| {sec['title']} | `results/{sec['name']}/` | {', '.join('`' + s.split(' ')[0] + '`' for s in sec['scripts'])} | {n} |")

    script_rows = '\n'.join(f'| `{s}` | {w} |' for s, w in SCRIPTS)
    with open(os.path.join(DST, 'README.md'), 'w', encoding='utf-8') as fh:
        fh.write(f"""# netGS in maize: code and result files

This folder holds the code and the result files behind the report
"netGS in maize: what we did, where each number comes from, and why it comes out that way".
`results/` follows the report's sections; each sub-folder has a README that names the columns behind the numbers in the report.
`code/` holds every script with its curated inputs and runs as it is.
`results/` is rebuilt from `code/output/` by `python make_package.py` (run in `code/`).

## Report section to files

| Report section | Result files | Scripts | Files |
|---|---|---|---|
{chr(10).join(rows)}

## Scripts

| Script | What it does |
|---|---|
{script_rows}

## Editable inputs (`code/curation/`)

| File | Content |
|---|---|
| `metabolite_mapping.csv` | Zhang metabolite name to C4GEM name |
| `anchors.csv` | anchor per metabolite (millimoles per gram dry weight) and its source |
| `bof_tableIII.csv` | C4GEM biomass reaction from Dal'Molin et al. 2010, Table III |
| `molecular_weights.csv` | molecular weights of the 24 new metabolites |
| `gapfill_reactions.json` | reactions added to the model |
| `broken_reactions.csv` | reactions switched off |

{RUN}""")
    print('written', res, 'and README.md')


if __name__ == '__main__':
    main()
