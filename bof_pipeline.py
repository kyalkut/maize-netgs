import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
"""
Line- and condition-specific biomass objective functions (BOFs) for maize
from relative metabolomics (Zhang et al. 2021) and C4GEM (Dal'Molin et al. 2010).

Pipeline (each step is a function below; run all with `python bof_pipeline.py`):

  1  load      Read C4GEM (SBML), Zhang Table S2/S3 (metabolite info + per-line levels), Table III BOF
  2  map       Apply the curated metabolite mapping (curation/metabolite_mapping.csv)
  3  produce   Producibility test: can the model synthesise each mapped metabolite from CO2, light, N, S, P?
  4  clean     Select features, collapse duplicates, cap outliers, impute missing values
  5  relative  Relative factors r_ij = 2^(x_ij - mean_WW(x_j))   (common WW anchor for both conditions)
  6  anchor    Absolute free-pool anchors per metabolite (curation/anchors.csv; placeholders flagged)
  7  twopool   b_ij = b_bound_j + b_free_j * r_ij
  8  mass      Rescale fixed components so every line has the same total BOF mass as the reference BOF
  9  export    biomreaction_WW.csv / biomreaction_DS.csv (netGS format: metabolites x lines) + long table
  10 qc        Feasibility (max biomass per line with its BOF) and summary report

Only numpy / scipy / openpyxl are needed.
"""
import csv, html, json, os, pickle, re, sys
import numpy as np
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
P = dict(
    model=os.path.join(HERE, 'input', 'Maize_C4GEM_vs1.0.xml'),
    zhang=os.path.join(HERE, 'input', 'Zhang2021_AdditionalFile2.xlsx'),
    mapping=os.path.join(HERE, 'curation', 'metabolite_mapping.csv'),
    bof=os.path.join(HERE, 'curation', 'bof_tableIII.csv'),
    anchors=os.path.join(HERE, 'curation', 'anchors.csv'),
    gapfill=os.path.join(HERE, 'curation', 'gapfill_reactions.json'),
    broken=os.path.join(HERE, 'curation', 'broken_reactions.csv'),
    out=os.path.join(HERE, 'output'),
    cache=os.path.join(HERE, 'output', 'cache_zhang.pkl'),
)
# ---- parameters (the decisions to discuss) ---------------------------------------------
PARAM = dict(
    min_coverage=0.90,        # feature must be measured in >=90% of lines in BOTH conditions
    max_sd_log2=2.5,          # drop features whose spread across lines looks like noise
    cap_quantiles=(0.01, 0.99),  # winsorise log2 values per feature and condition
    default_free_frac_aa=0.05,   # PLACEHOLDER: free pool = 5% of Table III coefficient (amino acids)
    default_free_other=0.001,    # PLACEHOLDER: 1 umol/gDW for metabolites without literature value
    placeholder_scale=1.0,       # sensitivity: multiply all placeholder anchors
    reference_line='B73',
    light_uptake=100.0,       # arbitrary cap for LP tests (results are relative)
    atp_maintenance=30.0,     # Table III: 30 ATP -> ADP + Pi per unit biomass
)

# =========================================================================================
# 1 LOAD
# =========================================================================================
def read_sbml(path):
    t = open(path, encoding='utf-8').read()
    sp = re.findall(r'<species id="([^"]+)" name="([^"]+)" compartment="C_([^"]+)"([^/]*)/>', t)
    mets = [s[0] for s in sp]
    base = [re.sub(r'_(c|m|p|v|x|ext|acc|pmd|total)$', '', html.unescape(s[1])) for s in sp]
    comp = [s[2] for s in sp]
    boundary = np.array(['boundaryCondition="true"' in s[3] for s in sp])
    idx = {m: i for i, m in enumerate(mets)}
    rx = re.findall(r'<reaction id="([^"]+)"([^>]*)>(.*?)</reaction>', t, re.S)
    S = np.zeros((len(mets), len(rx))); lb = np.zeros(len(rx)); ub = np.full(len(rx), 1000.0); rids = []
    for j, (rid, attr, body) in enumerate(rx):
        rids.append(rid)
        lb[j] = 0.0 if 'reversible="false"' in attr else -1000.0
        for tag, sgn in (('listOfReactants', -1), ('listOfProducts', 1)):
            m = re.search(rf'<{tag}>(.*?)</{tag}>', body, re.S)
            if not m: continue
            for ref in re.findall(r'<speciesReference ([^>]*?)/?>', m.group(1)):
                s = re.search(r'species="([^"]+)"', ref).group(1)
                st = re.search(r'stoichiometry="([^"]+)"', ref)
                S[idx[s], j] += sgn * (float(st.group(1)) if st else 1.0)
    return dict(S=S, lb=lb, ub=ub, rxns=rids, mets=mets, base=base, comp=comp, boundary=boundary)

def _apply_gapfill_unused(model, path):
    """Add reactions missing from the published SBML (curation/gapfill_reactions.json):
    43 parsable reactions from Dal'Molin 2010 Table S1 + 2 manual links (H2O2 c->x, UMP kinase).
    Added as reversible; direction to be checked."""
    gf = json.load(open(path)); idx = {s: i for i, s in enumerate(model['mets'])}; cols = []
    for r in gf:
        c = np.zeros(len(model['mets']))
        for co, s in r['L']: c[idx[s]] -= co
        for co, s in r['R']: c[idx[s]] += co
        cols.append(c)
    model['S'] = np.hstack([model['S'], np.array(cols).T])
    model['lb'] = np.append(model['lb'], [-1000.0] * len(cols)); model['ub'] = np.append(model['ub'], [1000.0] * len(cols))
    model['rxns'] = model['rxns'] + ['GF_' + r['id'] for r in gf]
    return len(gf)

def species_for(model, name, prefer=('c', 'p', 'm', 'x', 'v')):
    cands = {model['comp'][i]: i for i, b in enumerate(model['base']) if b == name}
    for c in prefer:
        if c in cands: return cands[c]
    return None

def load_zhang(path, cache):
    if os.path.exists(cache):
        return pickle.load(open(cache, 'rb'))
    import openpyxl
    print('  reading Zhang xlsx (one-time, a few minutes)...')
    wb = openpyxl.load_workbook(path, read_only=True)
    s2 = {}
    for r in wb['Table S2'].iter_rows(min_row=8, values_only=True):
        if r[0]: s2[r[0]] = dict(name=r[10], cls=r[11], log2fc=r[5], fdr=r[7], mass=r[8])
    s1 = {}
    for r in wb['Table S1'].iter_rows(min_row=7, values_only=True):
        if r[0] is not None: s1[str(r[0])] = dict(dry=r[4] if isinstance(r[4], (int, float)) else None,
                                                   fresh=r[5] if isinstance(r[5], (int, float)) else None)
    it = wb['Table S3'].iter_rows(values_only=True)
    for _ in range(3): next(it)
    lines = [str(x) for x in next(it)[3:]]
    levels, dr = {}, {}
    for r in it:
        if not r[0]: continue
        levels[(r[0], r[1])] = np.array([float(v) if isinstance(v, (int, float)) else np.nan for v in r[3:]])
        dr[r[0]] = (r[2] == 'Yes')
    data = dict(s2=s2, s1=s1, lines=lines, levels=levels, dr=dr)
    pickle.dump(data, open(cache, 'wb'))
    return data

def read_csv(path):
    return list(csv.DictReader(open(path, encoding='utf-8')))

# =========================================================================================
# 3 PRODUCIBILITY
# =========================================================================================
def autotrophic_bounds(model):
    """Uptake: CO2, light, nitrate, NH3, H2S, sulfate, phosphate, water, O2. No sugar uptake."""
    S = model['S'].copy(); lb = model['lb'].copy(); ub = model['ub'].copy(); rx = list(model['rxns'])
    keep = ~model['boundary']                     # boundary species are not mass-balanced
    # sulfate_ext is not flagged boundary in the file -> release it; add phosphate uptake (missing in file)
    i_sulf = model['mets'].index('S_Sulfate_ext'); keep[i_sulf] = False
    i_pi = species_for(model, 'Orthophosphate', ('c',))
    col = np.zeros((S.shape[0], 1)); col[i_pi] = 1; S = np.hstack([S, col]); lb = np.append(lb, 0); ub = np.append(ub, 1000); rx.append('EX_Pi_added')
    # block sugar uptake (Glucose_ext -> c); sucrose export only
    j = rx.index('R_Ex8'); ub[j] = 0
    # light cap
    j = rx.index('R_hv_total'); ub[j] = PARAM['light_uptake']
    return dict(S=S[keep], lb=lb, ub=ub, rxns=rx, keep=keep)

def max_production(net, model, met_index, coef_vector=None):
    """Maximise a drain on one metabolite (or a full BOF if coef_vector given)."""
    S, lb, ub = net['S'], net['lb'], net['ub']
    drain = np.zeros(model['S'].shape[0])
    if coef_vector is None: drain[met_index] = -1
    else: drain = coef_vector
    d = drain[net['keep']].reshape(-1, 1)
    A = np.hstack([S, d]); c = np.zeros(A.shape[1]); c[-1] = -1
    bounds = list(zip(lb, ub)) + [(0, 1000)]
    res = linprog(c, A_eq=A, b_eq=np.zeros(A.shape[0]), bounds=bounds, method='highs')
    return -res.fun if res.status == 0 else 0.0

# =========================================================================================
# MAIN
# =========================================================================================
def main():
    # optional overrides for sensitivity runs
    if os.environ.get('BOF_OUT'): P['out'] = os.environ['BOF_OUT']
    if os.environ.get('BOF_PH_SCALE'): PARAM['placeholder_scale'] = float(os.environ['BOF_PH_SCALE'])
    os.makedirs(P['out'], exist_ok=True)
    log = []
    def say(*a):
        s = ' '.join(str(x) for x in a); print(s); log.append(s)

    say('STEP 1  load')
    model = read_sbml(P['model'])
    from model_curation import curate
    cur = curate(model, P['model'], P['gapfill'], P['broken'])
    say(f'  curation: +{cur["gapfill"]} gap-fill reactions (Table S1 + manual), {cur["broken"]} broken reactions blocked, '
        f'{len(cur["direction_changes"])} reactions set irreversible by EC rules')
    with open(os.path.join(P['out'], 'model_curation_log.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['reaction', 'change'])
        for r in json.load(open(P['gapfill'])): w.writerow(['GF_' + r['id'], 'added: ' + r.get('source', '') + ('' if r.get('reversible') else ' (irreversible)')])
        for r in read_csv(P['broken']): w.writerow([r['reaction'], 'blocked: ' + r['problem']])
        for r, why in cur['direction_changes']: w.writerow([r, 'direction: ' + why])
    Z = load_zhang(P['zhang'], P['cache'])
    bof = read_csv(P['bof'])
    lines = Z['lines']; nL = len(lines)
    say(f'  model: {model["S"].shape[0]} species x {model["S"].shape[1]} reactions; Zhang: {len(Z["s2"])} features, {nL} lines; Table III: {len(bof)} components')
    base_mass = sum(float(b['coef_mmol_per_gDW']) * float(b['mw_g_per_mol']) for b in bof)
    say(f'  Table III BOF mass = {base_mass:.1f} mg per unit biomass (not 1000 mg: Table III covers only part of dry mass)')

    say('STEP 2  map')
    mapping = [m for m in read_csv(P['mapping']) if m['include'].strip().lower() == 'yes']
    say(f'  {len(mapping)} curated Zhang features mapped to {len({m["model_metabolite"] for m in mapping})} model metabolites')

    say('STEP 3  producibility (autotrophic: CO2 + light + N/S/P)')
    net = autotrophic_bounds(model)
    bof_ok = {}
    for b in bof:
        i = species_for(model, b['model_metabolite'], (b['compartment'], 'c', 'p'))
        bof_ok[b['component']] = (i is not None) and max_production(net, model, i) > 1e-6
    say(f'  Table III components producible: {sum(bof_ok.values())}/{len(bof)}; not producible: {[k for k,v in bof_ok.items() if not v]}')
    prod = {}
    for mm in sorted({m['model_metabolite'] for m in mapping}):
        i = species_for(model, mm)
        prod[mm] = max_production(net, model, i) if i is not None else 0.0
    producible = {k for k, v in prod.items() if v > 1e-6}
    say(f'  mapped metabolites producible: {len(producible)}/{len(prod)}; not producible: {sorted(set(prod)-producible)}')

    say('STEP 4  clean data')
    tab3 = {b['model_metabolite']: b for b in bof}
    chosen, dropped = {}, []
    for m in mapping:
        fid, mm = m['zhang_id'], m['model_metabolite']
        if mm not in producible: dropped.append((fid, mm, 'not producible')); continue
        ww, ds = Z['levels'].get((fid, 'WW')), Z['levels'].get((fid, 'DS'))
        if ww is None or ds is None: dropped.append((fid, mm, 'not measured in both conditions')); continue
        cov = min(np.isfinite(ww).mean(), np.isfinite(ds).mean())
        sd = max(np.nanstd(ww), np.nanstd(ds))
        if cov < PARAM['min_coverage']: dropped.append((fid, mm, f'coverage {cov:.2f}')); continue
        if sd > PARAM['max_sd_log2']: dropped.append((fid, mm, f'noisy (SD log2 {sd:.2f})')); continue
        score = (cov, np.nanmean(np.r_[ww, ds]))   # duplicates: most complete, then most intense
        if mm not in chosen or score > chosen[mm]['score']:
            if mm in chosen: dropped.append((chosen[mm]['fid'], mm, 'duplicate feature (lower rank)'))
            chosen[mm] = dict(fid=fid, score=score, ww=ww.copy(), ds=ds.copy(), name=m['zhang_name'],
                              mass=float(m['exact_mass']) if m['exact_mass'] else None)
        else:
            dropped.append((fid, mm, 'duplicate feature (lower rank)'))
    lo, hi = PARAM['cap_quantiles']
    for mm, c in chosen.items():
        for k in ('ww', 'ds'):
            v = c[k]; ql, qh = np.nanquantile(v, lo), np.nanquantile(v, hi)
            v = np.clip(v, ql, qh); v[~np.isfinite(v)] = np.nanmean(v); c[k] = v
    say(f'  kept {len(chosen)} metabolites; dropped {len(dropped)} features')
    for d in dropped: say(f'    - {d[0]:16} {d[1]:35} {d[2]}')

    say('STEP 5  relative factors (anchor = WW mean, same for both conditions)')
    rel = {}
    for mm, c in chosen.items():
        mu = np.mean(c['ww'])
        rel[mm] = dict(WW=2 ** (c['ww'] - mu), DS=2 ** (c['ds'] - mu))
    for mm in ('L-Tryptophan', 'L-Methionine', 'L-Glutamate'):
        if mm in rel: say(f'  {mm:16} median r  WW={np.median(rel[mm]["WW"]):.2f}  DS={np.median(rel[mm]["DS"]):.2f}')

    say('STEP 6  anchors (absolute free-pool amounts, mmol/gDW)')
    anchors = {a['model_metabolite']: a for a in read_csv(P['anchors'])} if os.path.exists(P['anchors']) else {}
    free, src = {}, {}
    for mm in chosen:
        a = anchors.get(mm)
        if a and a.get('free_mmol_per_gDW', '').strip():
            free[mm] = float(a['free_mmol_per_gDW']); src[mm] = a.get('source', 'curated')
        elif mm in tab3:
            free[mm] = PARAM['placeholder_scale'] * PARAM['default_free_frac_aa'] * float(tab3[mm]['coef_mmol_per_gDW']); src[mm] = 'PLACEHOLDER (5% of Table III)'
        else:
            free[mm] = PARAM['placeholder_scale'] * PARAM['default_free_other']; src[mm] = 'PLACEHOLDER (1 umol/gDW)'
    n_ph = sum(s.startswith('PLACEHOLDER') for s in src.values())
    say(f'  {len(free)} anchors; {n_ph} are placeholders -> replace with literature values in curation/anchors.csv')
    # write/refresh anchors template
    if not os.environ.get('BOF_OUT'):
      with open(P['anchors'], 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow(['model_metabolite', 'zhang_name', 'in_tableIII', 'free_mmol_per_gDW', 'source'])
        for mm in sorted(chosen):
            a = anchors.get(mm, {})
            w.writerow([mm, chosen[mm]['name'], 'yes' if mm in tab3 else 'no',
                        a.get('free_mmol_per_gDW', '') if a else '', a.get('source', '') if a else ''])

    say('STEP 7  two-pool coefficients  b = b_bound + b_free * r')
    comps = [b['model_metabolite'] for b in bof] + sorted(m for m in chosen if m not in tab3)
    mw = {b['model_metabolite']: float(b['mw_g_per_mol']) for b in bof}
    for mm in chosen:
        if mm not in mw: mw[mm] = chosen[mm]['mass'] or 150.0
    bound = np.array([float(tab3[c]['coef_mmol_per_gDW']) if c in tab3 else 0.0 for c in comps])
    is_var = np.array([c in chosen for c in comps])
    out = {}
    for cond in ('WW', 'DS'):
        B = np.tile(bound.reshape(-1, 1), (1, nL))
        for k, c in enumerate(comps):
            if c in chosen: B[k] = bound[k] + free[c] * rel[c][cond]
        out[cond] = B

    say('STEP 8  mass renormalisation (fixed part rescaled, measured differences kept)')
    MW = np.array([mw[c] for c in comps])
    ref_B = bound + np.array([free.get(c, 0.0) if c in chosen else 0.0 for c in comps])   # average line, r = 1
    M0 = float(ref_B @ MW)
    fixed_scale = {}
    for cond in ('WW', 'DS'):
        B = out[cond]
        var_mass = (B[is_var] - bound[is_var].reshape(-1, 1)).T @ MW[is_var]       # free-pool mass per line
        fixed_mass = bound @ MW
        s = (M0 - var_mass) / fixed_mass                                             # scale for fixed part
        fixed_scale[cond] = s
        B = bound.reshape(-1, 1) * s + np.where(is_var.reshape(-1, 1), B - bound.reshape(-1, 1), 0)
        out[cond] = B
        mass = B.T @ MW
        say(f'  {cond}: target mass {M0:.1f} mg; per-line mass range {mass.min():.2f}-{mass.max():.2f}; fixed-part scale {s.min():.4f}-{s.max():.4f}')

    say('STEP 9  export')
    for cond in ('WW', 'DS'):
        with open(os.path.join(P['out'], f'biomreaction_{cond}.csv'), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['model_metabolite', 'compartment', 'varies_by_line'] + lines)
            for k, c in enumerate(comps):
                comp = tab3[c]['compartment'] if c in tab3 else model['comp'][species_for(model, c)]
                w.writerow([c, comp, 'yes' if is_var[k] else 'no'] + [f'{x:.6g}' for x in out[cond][k]])
    with open(os.path.join(P['out'], 'bof_components.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['model_metabolite', 'in_tableIII', 'varies', 'zhang_feature', 'bound_mmol_per_gDW', 'free_anchor_mmol_per_gDW', 'anchor_source', 'mw'])
        for k, c in enumerate(comps):
            w.writerow([c, 'yes' if c in tab3 else 'no', 'yes' if is_var[k] else 'no', chosen[c]['fid'] if c in chosen else '',
                        bound[k], free.get(c, '') if c in chosen else '', src.get(c, '') if c in chosen else '', MW[k]])
    say(f'  {len(comps)} BOF components ({int(is_var.sum())} vary by line) x {nL} lines, WW and DS')

    # ---- rescaling trace: every step from log2 intensity to final coefficient -------------------
    var_comps = [c for c in comps if c in chosen]
    for cond in ('WW', 'DS'):
        for name, getter in (('relative_factors', lambda c: rel[c][cond]),
                             ('free_pool_mmol', lambda c: free[c] * rel[c][cond])):
            with open(os.path.join(P['out'], f'{name}_{cond}.csv'), 'w', newline='') as f:
                w = csv.writer(f); w.writerow(['model_metabolite', 'zhang_feature'] + lines)
                for c in var_comps:
                    w.writerow([c, chosen[c]['fid']] + [f'{x:.6g}' for x in getter(c)])
    with open(os.path.join(P['out'], 'rescaling_trace.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['condition', 'line', 'model_metabolite', 'zhang_feature',
                    '1_log2_intensity', '2_mean_log2_WW', '3_relative_factor', '4_free_anchor_mmol',
                    '5_free_part_mmol', '6_bound_part_tableIII', '7_fixed_part_scale', '8_bound_part_scaled', '9_final_coefficient'])
        for cond in ('WW', 'DS'):
            for c in var_comps:
                k = comps.index(c); x = chosen[c]['ww' if cond == 'WW' else 'ds']; mu = float(np.mean(chosen[c]['ww']))
                for li, line in enumerate(lines):
                    r = rel[c][cond][li]; fp = free[c] * r; sc = fixed_scale[cond][li]
                    w.writerow([cond, line, c, chosen[c]['fid'], f'{x[li]:.4f}', f'{mu:.4f}', f'{r:.5g}', f'{free[c]:.5g}',
                                f'{fp:.5g}', f'{bound[k]:.5g}', f'{sc:.5f}', f'{bound[k] * sc:.5g}', f'{out[cond][k, li]:.6g}'])
    say('  rescaling trace: rescaling_trace.csv (long), relative_factors_*.csv, free_pool_mmol_*.csv (metabolites x lines)')

    say('STEP 10 QC: feasibility of every line-specific BOF (single-cell network)')
    feas = {}
    i_atp, i_adp, i_pi = (species_for(model, n, ('c',)) for n in ('ATP', 'ADP', 'Orthophosphate'))
    for cond in ('WW', 'DS'):
        vals = []
        for li in range(nL):
            v = np.zeros(model['S'].shape[0])
            for k, c in enumerate(comps):
                i = species_for(model, c, (tab3[c]['compartment'], 'c', 'p') if c in tab3 else ('c', 'p', 'm', 'x', 'v'))
                v[i] -= out[cond][k, li]
            v[i_atp] -= PARAM['atp_maintenance']; v[i_adp] += PARAM['atp_maintenance']; v[i_pi] += PARAM['atp_maintenance']
            vals.append(max_production(net, model, None, v))
        feas[cond] = np.array(vals)
        f = feas[cond]
        say(f'  {cond}: max biomass flux per line: min {f.min():.4g}, mean {f.mean():.4g}, max {f.max():.4g}; infeasible lines: {(f<=1e-9).sum()}')
    ref = lines.index(PARAM['reference_line'])
    say(f'  reference {PARAM["reference_line"]}: WW {feas["WW"][ref]:.4g}, DS {feas["DS"][ref]:.4g}')
    # energy-loop check: can ATP be regenerated with zero light?
    saved = PARAM['light_uptake']; PARAM['light_uptake'] = 0.0
    net0 = autotrophic_bounds(model); a = np.zeros(model['S'].shape[0]); a[i_atp] = -1; a[i_adp] = 1; a[i_pi] = 1
    free_atp = max_production(net0, model, None, a); PARAM['light_uptake'] = saved
    say(f'  energy-loop check: ATP regenerated without light = {free_atp:.0f}  ->  '
        + ('LOOPS PRESENT: feasibility above = connectivity only; absolute biomass values (s_model) not yet meaningful'
           if free_atp > 1e-6 else 'no loops'))

    # figure: per-line relative change of each varying coefficient (vs reference BOF), WW vs DS
    try:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        var_idx = [k for k in range(len(comps)) if is_var[k]]
        order = sorted(var_idx, key=lambda k: np.median(out['DS'][k] / ref_B[k]))
        fig, ax = plt.subplots(figsize=(8, max(4, 0.28 * len(order))))
        pos = np.arange(len(order))
        for cond, off, col in (('WW', -0.18, '#2a78d6'), ('DS', 0.18, '#eb6834')):
            data = [np.log2(out[cond][k] / ref_B[k]) for k in order]
            bp_ = ax.boxplot(data, positions=pos + off, widths=0.32, orientation='horizontal', patch_artist=True, showfliers=False)
            for b in bp_['boxes']: b.set_facecolor(col); b.set_alpha(0.7)
        ax.axvline(0, color='grey', lw=0.8)
        ax.set_yticks(pos); ax.set_yticklabels([comps[k] for k in order], fontsize=8)
        ax.set_xlabel('log2( line coefficient / reference coefficient )   [blue = WW, orange = DS]')
        ax.set_title('Line-specific BOF coefficients (two-pool, placeholder anchors)', fontsize=10)
        fig.tight_layout(); fig.savefig(os.path.join(P['out'], 'qc_bof_coefficients.png'), dpi=150)
        say('  figure: output/qc_bof_coefficients.png')
    except Exception as e:
        say(f'  (figure skipped: {e})')
    with open(os.path.join(P['out'], 'feasibility_per_line.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['line', 'max_biomass_WW', 'max_biomass_DS', 'dry_mass_DS_g'])
        for li, l in enumerate(lines):
            w.writerow([l, f'{feas["WW"][li]:.6g}', f'{feas["DS"][li]:.6g}', Z['s1'].get(l, {}).get('dry', '')])
    open(os.path.join(P['out'], 'run_log.txt'), 'w').write('\n'.join(log))
    json.dump(dict(param=PARAM, n_components=len(comps), n_varying=int(is_var.sum()), M0=M0,
                   kept=sorted(chosen), dropped=dropped, producible=sorted(producible),
                   not_producible=sorted(set(prod) - producible), tableIII_not_producible=[k for k, v in bof_ok.items() if not v]),
              open(os.path.join(P['out'], 'run_summary.json'), 'w'), indent=1, default=str)

if __name__ == '__main__':
    main()
