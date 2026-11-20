"""Compare the two BOF variants (both normalised to 1 g DW):
  with24 = output/        Table III + 7 amino acids (Shu anchors) + 24 metabolites with placeholder anchors (61 components)
  no24   = output/no24/   Table III + 7 amino acids only (37 components; the 24 left out, BOF_NO_PH=1)
Part 1 (BOF): coefficients of the shared components, line max biomass, spread, ranking, drought shift.
Part 2 (netGS step 1, if the flux files exist): per-reaction correlation across lines, variation beyond growth, flux levels.
Writes output/variant_comparison.json and output/variant_comparison_reactions_<tag>.csv."""
import os, csv, json, numpy as np
from scipy.stats import spearmanr, pearsonr
HERE = os.path.dirname(os.path.abspath(__file__))
D = {'with24': os.path.join(HERE, 'output'), 'no24': os.path.join(HERE, 'output', 'no24')}

def bof(var, cond):
    rows = list(csv.reader(open(os.path.join(D[var], f'biomreaction_{cond}.csv'))))
    return rows[0][3:], {r[0]: np.array(r[3:], float) for r in rows[1:]}, {r[0]: r[2] for r in rows[1:]}

def feas(var):
    rows = list(csv.DictReader(open(os.path.join(D[var], 'feasibility_per_line.csv'))))
    return [r['line'] for r in rows], {c: np.array([float(r[f'max_biomass_{c}']) for r in rows]) for c in ('WW', 'DS')}

out = {}
# ---------------- Part 1: BOF ----------------
L, F1 = feas('with24'); L2, F0 = feas('no24'); assert L == L2
p1 = {}
for c in ('WW', 'DS'):
    a, b = F1[c], F0[c]
    p1[c] = dict(max_biomass_with24=[float(a.min()), float(a.max())], max_biomass_no24=[float(b.min()), float(b.max())],
                 spread_pct_with24=float(100 * (a.max() - a.min()) / a.mean()), spread_pct_no24=float(100 * (b.max() - b.min()) / b.mean()),
                 spearman_line_ranking=float(spearmanr(a, b)[0]), B73_with24=float(a[L.index('B73')]), B73_no24=float(b[L.index('B73')]))
# drought shift in max biomass per line (DS - WW), and its agreement between variants
d1, d0 = F1['DS'] - F1['WW'], F0['DS'] - F0['WW']
p1['drought_shift'] = dict(median_with24=float(np.median(d1)), median_no24=float(np.median(d0)),
                           sd_with24=float(np.std(d1)), sd_no24=float(np.std(d0)), spearman=float(spearmanr(d1, d0)[0]))
# shared components: the 37 Table III ones (7 vary by line in both)
coef = {}
for c in ('WW', 'DS'):
    lines, B1, V1 = bof('with24', c); _, B0, V0 = bof('no24', c)
    rows = []
    for k in B0:
        r_ = spearmanr(B1[k], B0[k])[0] if V0[k] == 'yes' else np.nan
        rows.append(dict(component=k, varies=V0[k], median_with24=float(np.median(B1[k])), median_no24=float(np.median(B0[k])),
                         ratio_no24_over_with24=float(np.median(B0[k] / B1[k])), spearman_lines=float(r_)))
    coef[c] = rows
    if c == 'DS':
        with open(os.path.join(D['with24'], 'variant_comparison_bof_DS.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
            for r in rows: w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
rat = np.array([r['ratio_no24_over_with24'] for r in coef['DS']])
p1['shared_37_DS'] = dict(ratio_no24_over_with24_range=[float(rat.min()), float(rat.max())],
                          spearman_lines_varying_7=[r['spearman_lines'] for r in coef['DS'] if r['varies'] == 'yes'])
mass24 = {}
for var in D:
    comp = list(csv.DictReader(open(os.path.join(D[var], 'bof_components.csv'))))
    lines, B, V = bof(var, 'DS'); mw = {c['model_metabolite']: float(c['mw']) for c in comp}
    i = lines.index('B73')
    tab = sum(B[c['model_metabolite']][i] * mw[c['model_metabolite']] for c in comp if c['in_tableIII'] == 'yes')
    mass24[var] = dict(n_components=len(comp), n_varying=sum(c['varies'] == 'yes' for c in comp),
                       B73_DS_mass_tableIII_part_mg=float(tab), B73_DS_mass_other_mg=float(1000 - tab))
p1['composition'] = mass24
out['bof'] = p1

# ---------------- Part 2: netGS step 1 ----------------
def flux(var, tag):
    p = os.path.join(D[var], f'netgs_fluxes_DS_{tag}.csv')
    if not os.path.exists(p): return None
    rows = list(csv.reader(open(p)))
    return rows[0][2:], {r[0]: np.array(r[2:], float) for r in rows[1:]}, {r[0]: float(r[1]) for r in rows[1:]}

def beyond(Fd):
    bio = Fd['BIOMASS']; cv = {}
    for k, v in Fd.items():
        u = v / bio; m = np.nanmean(u)
        cv[k] = float(np.nanstd(u) / abs(m)) if abs(m) > 0 else np.nan
    return cv

p2 = {}
for tag in ('paper_cvxpy_ser', 'code_cvxpy_ser', 'paper_cvxpy_co', 'code_cvxpy_co'):
    a, b = flux('with24', tag), flux('no24', tag)
    if a is None or b is None: continue
    L1, F1_, R1 = a; L0, F0_, R0 = b
    if L1 != L0: continue   # run not finished yet
    shared = [k for k in F1_ if k in F0_]
    cv1, cv0 = beyond(F1_), beyond(F0_)
    rows = []
    for k in shared:
        x, y = F1_[k], F0_[k]
        r_ = pearsonr(x, y)[0] if np.std(x) > 0 and np.std(y) > 0 else np.nan
        rows.append(dict(reaction=k, v_ref_with24=R1[k], v_ref_no24=R0[k], r_lines=float(r_),
                         level_ratio_no24_over_with24=float(np.median(y / x)) if np.all(x != 0) else np.nan,
                         cv_beyond_growth_with24=cv1[k], cv_beyond_growth_no24=cv0[k]))
    with open(os.path.join(D['with24'], f'variant_comparison_reactions_{tag}.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
        for r in rows: w.writerow({k: (f'{v:.6g}' if isinstance(v, float) else v) for k, v in r.items()})
    rr = np.array([r['r_lines'] for r in rows]); lr = np.array([r['level_ratio_no24_over_with24'] for r in rows])
    n1 = sum(v > 0.05 for v in cv1.values() if np.isfinite(v)); n0 = sum(v > 0.05 for v in cv0.values() if np.isfinite(v))
    only1 = sorted(k for k in F1_ if k not in F0_)
    p2[tag] = dict(n_active_with24=len(F1_), n_active_no24=len(F0_), n_shared=len(shared), n_only_with24=len(only1),
                   only_no24=sorted(k for k in F0_ if k not in F1_),
                   r_lines_median=float(np.nanmedian(rr)), r_lines_min=float(np.nanmin(rr)), n_r_below_0_9=int(np.sum(rr < 0.9)),
                   level_ratio_median=float(np.nanmedian(lr)), level_ratio_range=[float(np.nanmin(lr)), float(np.nanmax(lr))],
                   n_reactions_cv_gt_5pct_with24=int(n1), n_reactions_cv_gt_5pct_no24=int(n0),
                   median_cv_beyond_growth_with24=float(np.nanmedian(list(cv1.values()))), median_cv_beyond_growth_no24=float(np.nanmedian(list(cv0.values()))),
                   biomass_r_lines=float(pearsonr(F1_['BIOMASS'], F0_['BIOMASS'])[0]))
out['netgs'] = p2
json.dump(out, open(os.path.join(D['with24'], 'variant_comparison.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
