"""Copy every step's inputs and outputs into steps/<NN_step>/inputs and /outputs, with a README per step.
The pipeline itself still reads input/, curation/ and writes output/; steps/ is a browsable copy.
Rerun after any pipeline run:  python make_steps.py"""
import os, shutil

HERE = os.path.dirname(os.path.abspath(__file__)); DST = os.path.join(HERE, 'steps')
ZHANG = 'input/Zhang2021_AdditionalFile2.xlsx'

STEPS = [
    dict(name='01_model_curation',
         what="Repair C4GEM: +46 gap-fill reactions (43 from Dal'Molin 2010 Table S1), 12 broken reactions blocked, 28 directions set by EC rules.",
         scripts='model_curation.py (called by bof_pipeline.py, step 1)',
         inputs=['input/Maize_C4GEM_vs1.0.xml', 'curation/gapfill_reactions.json', 'curation/broken_reactions.csv'],
         outputs=['output/model_curation_log.csv']),
    dict(name='02_metabolite_selection',
         what='351 annotated Zhang features -> 52 name matches (50 metabolites) -> 35 producible -> 31 with clean data.',
         scripts='bof_pipeline.py, steps 1-4',
         inputs=[ZHANG, 'curation/metabolite_mapping.csv', 'curation/bof_tableIII.csv'],
         outputs=['output/run_log.txt', 'output/run_summary.json', 'output/bof_components.csv']),
    dict(name='03_rescaling',
         what='log2 data -> mmol/gDW: relative factor 2^(x - WW mean), x anchor, + Table III part, x mass factor s (536.3 mg).',
         scripts='bof_pipeline.py, steps 5-9',
         inputs=['curation/anchors.csv', 'curation/bof_tableIII.csv', 'output/bof_components.csv'],
         outputs=['output/relative_factors_WW.csv', 'output/relative_factors_DS.csv',
                  'output/free_pool_mmol_WW.csv', 'output/free_pool_mmol_DS.csv',
                  'output/rescaling_trace.csv', 'output/biomreaction_WW.csv', 'output/biomreaction_DS.csv',
                  'output/qc_bof_coefficients.png'],
         note='Zhang log2 values: see 02_metabolite_selection/inputs.'),
    dict(name='04_bof_feasibility_and_sensitivity',
         what='Max biomass of all 770 BOFs (FBA); placeholder anchors x0.1 and x10.',
         scripts='bof_pipeline.py, step 10; sensitivity: BOF_PH_SCALE=0.1 / 10 with BOF_OUT=output/sensitivity_x0.1 / _x10',
         inputs=['output/biomreaction_WW.csv', 'output/biomreaction_DS.csv', 'input/Maize_C4GEM_vs1.0.xml'],
         outputs=['output/feasibility_per_line.csv', 'output/sensitivity_x0.1', 'output/sensitivity_x10']),
    dict(name='05_netgs_step1',
         what='Per DS line: s_model, growth rescaling from dry mass, B73 reference (pFBA), QP closest to reference.',
         scripts='run_netgs_maize.py (netgs_step1.py; solver cvxpy/Clarabel, qp_admm.py = old ADMM), merge_netgs.py (figure)',
         inputs=['output/biomreaction_DS.csv', 'output/feasibility_per_line.csv', 'input/Maize_C4GEM_vs1.0.xml'],
         outputs=['output/netgs_fluxes_DS_paper.csv', 'output/netgs_summary_DS_paper.json', 'output/fig_netgs_fluxes.png',
                  'output/netgs_fluxes_DS_code_cvxpy.csv', 'output/netgs_summary_DS_code_cvxpy.json',
                  'output/netgs_fluxes_DS_code_cvxpy_ser.csv', 'output/netgs_summary_DS_code_cvxpy_ser.json',
                  'output/netgs_fluxes_DS_paper_cvxpy_ser.csv', 'output/netgs_summary_DS_paper_cvxpy_ser.json'],
         note='Dry mass per line: column dry_mass_DS_g of feasibility_per_line.csv (from Zhang Table S1). '
              'The ADMM run (netgs_fluxes_DS_paper.csv) has 67 lines that are infeasible on the B73-active network (see the README). '
              'Files ending _ser add the phosphorylated serine pathway (R01513, R04173, R00582) to that network; all 371 lines are then feasible. '
              'code_cvxpy = objective and solver type of the original 1_FluxDist.m; paper_cvxpy = Eq. 2 of the paper (both interior point, Clarabel).'),
    dict(name='06_data_analyses',
         what='G x E consistency per BOF metabolite (WW vs DS line correlation); ridge CV prediction of DS dry mass.',
         scripts='analysis_data.py, figures.py',
         inputs=['output/bof_components.csv', 'output/run_summary.json', 'curation/metabolite_mapping.csv'],
         outputs=['output/analysis_GxE_bof_metabolites.csv', 'output/analysis_summary.json',
                  'output/fig_GxE_consistency.png', 'output/fig_prediction_drymass.png'],
         note='Zhang data: see 02_metabolite_selection/inputs.'),
    dict(name='08_fluxsum',
         what='Flux-sum phi = 1/2 sum|S_ij v_j| of each measured metabolite per DS line (from step 05 fluxes); Spearman r with the measured DS log2 level across lines.',
         scripts='fluxsum.py (NETGS_TAG picks the flux file, default paper_cvxpy_ser)',
         inputs=['output/netgs_fluxes_DS_paper_cvxpy_ser.csv', 'output/biomreaction_DS.csv', 'output/bof_components.csv',
                 'curation/metabolite_mapping.csv', 'input/Maize_C4GEM_vs1.0.xml'],
         outputs=['output/fluxsum_DS_paper_cvxpy_ser.csv', 'output/fluxsum_vs_level_DS_paper_cvxpy_ser.csv',
                  'output/fluxsum_summary_DS_paper_cvxpy_ser.json', 'output/fig_fluxsum.png'],
         note='Zhang log2 values: see 02_metabolite_selection/inputs. Main variant excludes the biomass reaction (circular); '
              'rho_level_phi_with_bof includes it.'),
    dict(name='07_validation_arabidopsis',
         what='Python port vs published netGS on the Arabidopsis example (s_model identical, 0.0035036); paper (Eq. 2) vs code objective: same flux patterns, r = 1.00 (cvxpy).',
         scripts='validation/validate_arabidopsis.py, validation/compare_objectives.py (NETGS_REPO = netGS repo netGS/ folder)',
         inputs=[],
         outputs=['validation/arab_fluxes_paper.npy', 'validation/arab_fluxes_code.npy'],
         note='Inputs (biomreaction.csv, biomass_optN.csv, model) come from github.com/Hao-Tong/netGS and are not stored here.'),
]

def copy(rel, sub):
    src = os.path.join(HERE, rel)
    if not os.path.exists(src): return None
    dst = os.path.join(sub, os.path.basename(rel))
    if os.path.isdir(src): shutil.copytree(src, dst)
    else: shutil.copy2(src, dst)
    return os.path.basename(rel) + ('/' if os.path.isdir(src) else '')

if os.path.isdir(DST): shutil.rmtree(DST)
index = ['# Pipeline steps', '', 'Copies of each step\'s inputs and outputs. Made by `make_steps.py`; the scripts still use input/, curation/ and output/.', '']
for s in STEPS:
    d = os.path.join(DST, s['name']); got = {}
    for kind in ('inputs', 'outputs'):
        sub = os.path.join(d, kind); os.makedirs(sub)
        got[kind] = [(rel, copy(rel, sub)) for rel in s[kind]]
    md = [f"# {s['name']}", '', s['what'], '', f"**Scripts:** {s['scripts']}", '']
    for kind in ('inputs', 'outputs'):
        md.append(f'**{kind.capitalize()}:**')
        md += [f'- `{kind}/{b}` (from `{rel}`)' if b else f'- `{rel}`: not found, run the step first' for rel, b in got[kind]] or ['- none stored here']
        md.append('')
    if s.get('note'): md += [s['note'], '']
    open(os.path.join(d, 'README.md'), 'w', encoding='utf-8').write('\n'.join(md))
    index.append(f"- `{s['name']}/`: {s['what']}")
open(os.path.join(DST, 'README.md'), 'w', encoding='utf-8').write('\n'.join(index) + '\n')
print('written', DST)
