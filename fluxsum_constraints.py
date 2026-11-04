"""Alternative to data-driven BOF components: metabolomics enters netGS step 1 through flux-sums.
  BOF: Table III (bound parts, mmol/gDW), the same for every line; lines differ by growth (dry mass) and flux-sums.
  Reference B73: pFBA with that BOF; active reactions kept, flux signs fixed to the reference (|S_ij v_j| becomes linear).
  Flux-sum of metabolite i (all compartments): phi_i = 1/2 sum_j |S_ij| s_j v_j,  s_j = sign(v_ref,j).
  Target per line: phi_i* = phi_i,ref * (bio_line / bio_ref) * 2^(x_line - x_B73)   (turnover per unit growth ~ level)
  QP per line (paper Eq. 2 in y = v / v_ref):  min 1/2 ||y - 1||^2 + lam/2 * sum_i (phi_i / phi_i* - 1)^2
      s.t. S v = 0, bounds, sign(v) = sign(v_ref), biomass in the growth window.
Cross-check (not circular): constrain one half of the metabolites, test the other half's phi against its measured level.
Env: FS_LAMBDA (default 10), FS_REF (fba / pfba), NETGS_N (lines)."""
import os, sys, csv, json, numpy as np, cvxpy as cp, scipy.sparse as sp
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
from model_curation import curate
from netgs_step1 import fba
from scipy.optimize import linprog
LAM = float(os.environ.get('FS_LAMBDA', '10')); NMAX = int(os.environ.get('NETGS_N', '100000')); OUT = bp.P['out']
m = bp.read_sbml(bp.P['model']); curate(m, bp.P['model'], bp.P['gapfill'], bp.P['broken'])
bp.PARAM['light_uptake'] = 100.0; net = bp.autotrophic_bounds(m)
base = np.array(m['base'])[net['keep']]
# Table III BOF (bound parts) + ATP maintenance
comp = {c['model_metabolite']: c for c in csv.DictReader(open(os.path.join(OUT, 'bof_components.csv')))}
rows = list(csv.reader(open(os.path.join(OUT, 'biomreaction_DS.csv'))))
bof = np.zeros(m['S'].shape[0])
for r in rows[1:]:
    b = float(comp[r[0]]['bound_mmol_per_gDW'] or 0)
    if b: bof[bp.species_for(m, r[0], (r[1], 'c', 'p', 'm', 'x', 'v'))] -= b
for nm, s in (('ATP', -1), ('ADP', 1), ('Orthophosphate', 1)): bof[bp.species_for(m, nm, ('c',))] += s * bp.PARAM['atp_maintenance']
bof = bof[net['keep']]
S = np.hstack([net['S'], bof[:, None]]); lb = np.append(net['lb'], 0); ub = np.append(net['ub'], 1000); jb = S.shape[1] - 1
c = np.zeros(S.shape[1]); c[jb] = 1
REF = os.environ.get("FS_REF", "pfba")  # reference flux: "pfba" (as run_netgs_maize.py) or "fba" (as the original 1_FluxDist.m)
vref, zmax = fba(S, lb, ub, c, pfba=(REF == "pfba")); nz = np.where(np.abs(vref) > 1e-5)[0]; vr = vref[nz]; kb = int(np.where(nz == jb)[0][0])
Sn = S[:, nz]; sgn = np.sign(vr)
# growth targets (same rule as netGS step 1; one BOF, so s_model = its max biomass)
Z = bp.load_zhang(bp.P['zhang'], bp.P['cache']); zl = Z['lines']
lines = rows[0][3:]; dm = np.array([Z['s1'].get(l, {}).get('dry') or np.nan for l in lines], float)
keep = np.where(np.isfinite(dm))[0][:NMAX]; L = [lines[k] for k in keep]; li = [zl.index(l) for l in L]
target = (zmax - 1.1e-4) * dm[keep] / np.nanmax(dm[keep]); half = 1.645 * np.std(target, ddof=1) / np.sqrt(len(target))
# measured metabolites with turnover in the reference: phi = g . v
mets, G, X = [], [], []
seen = set()
for mp in bp.read_csv(bp.P['mapping']):
    name, fid = mp['model_metabolite'], mp['zhang_id']
    if name in seen: continue
    rws = np.where(base == name)[0]; lev = Z['levels'].get((fid, 'DS'))
    if len(rws) == 0 or lev is None: continue
    g = 0.5 * np.abs(Sn[rws]).sum(0) * sgn; g[kb] = 0           # biomass column left out (no data in the BOF anyway)
    if g @ vr <= 1e-9: continue
    x = lev[li].astype(float); lo, hi = np.nanpercentile(x, [1, 99]); x = np.clip(x, lo, hi)
    if np.isfinite(x).mean() < 0.9: continue
    seen.add(name); mets.append(name); G.append(g); X.append(x)
G = np.array(G); X = np.array(X); phi_ref = G @ vr
xb = X[:, L.index('B73')] if 'B73' in L else np.nanmean(X, 1)
print(f'{len(nz)} active reactions; {len(mets)} measured metabolites with turnover in B73: {mets}')

def solve(i, use):
    """QP for line i; use = indices of metabolites whose flux-sum is constrained."""
    y = cp.Variable(len(nz)); v = cp.multiply(vr, y)
    cons = [sp.csr_matrix(Sn) @ v == 0, y >= 0, v[kb] >= target[i] - half, v[kb] <= target[i] + half]
    lo_, hi_ = lb[nz], ub[nz]
    cons += [v >= np.minimum(lo_, 0) * (sgn > 0) + lo_ * (sgn < 0), v <= hi_ * (sgn > 0) + np.maximum(hi_, 0) * (sgn < 0)]
    obj = 0.5 * cp.sum_squares(y - 1)
    u = [k for k in use if np.isfinite(X[k, i])]
    if u:
        tgt = phi_ref[u] * (target[i] / vr[kb]) * 2.0 ** (X[u, i] - xb[u])
        obj = obj + 0.5 * LAM * cp.sum_squares(cp.multiply(1.0 / tgt, G[u] @ v) - 1)
    pr = cp.Problem(cp.Minimize(obj), cons)
    try: pr.solve(solver=cp.CLARABEL)
    except cp.SolverError: pass
    ok = pr.status in ('optimal', 'optimal_inaccurate')
    return (np.asarray(v.value) if ok else np.full(len(nz), np.nan)), pr.status

def rho_tab(V, idx):
    out = {}
    for k in idx:
        phi = G[k] @ V; ok = np.isfinite(phi) & np.isfinite(X[k])
        out[mets[k]] = float(spearmanr(phi[ok], X[k, ok])[0]) if np.std(phi[ok]) > 1e-12 * np.mean(np.abs(phi[ok])) else np.nan
    return out

rng = np.random.default_rng(0); perm = rng.permutation(len(mets)); folds = [sorted(perm[::2]), sorted(perm[1::2])]
runs = {'none': [], 'all': list(range(len(mets))), 'foldA': folds[0], 'foldB': folds[1]}
V = {}; status = {}
for nm, use in runs.items():
    res = [solve(i, use) for i in range(len(keep))]
    V[nm] = np.array([r[0] for r in res]).T; status[nm] = {s: [r[1] for r in res].count(s) for s in set(r[1] for r in res)}
    print(nm, status[nm])
allm = range(len(mets))
tab = dict(none=rho_tab(V['none'], allm), all=rho_tab(V['all'], allm),
           held_out={**rho_tab(V['foldA'], folds[1]), **rho_tab(V['foldB'], folds[0])})
def med(d): v = np.array([x for x in d.values() if np.isfinite(x)]); return float(np.median(v)) if len(v) else None
# flux pattern beyond growth
def beyond(F):
    U = F / F[kb]; cv = np.nanstd(U, 1) / np.abs(np.nanmean(U, 1)); return float(np.nanmedian(cv)), int(np.sum(cv > 0.05))
summary = dict(lambda_=LAM, n_lines=len(keep), n_active=int(len(nz)), metabolites=mets, folds=[[mets[k] for k in f] for f in folds],
               status=status, median_rho=dict(none=med(tab['none']), all=med(tab['all']), held_out=med(tab['held_out'])),
               rho_per_metabolite=tab, beyond_growth=dict(none=beyond(V['none']), all=beyond(V['all'])),
               mean_abs_rel_dev_from_target_all=float(np.nanmedian(np.abs((G @ V['all']) / (phi_ref[:, None] * (target[None, :] / vr[kb]) * 2.0 ** (X - xb[:, None])) - 1))))
json.dump(summary, open(os.path.join(OUT, 'fluxsum_constraints_summary_DS.json'), 'w'), indent=1)
with open(os.path.join(OUT, 'fluxsum_constraints_fluxes_DS.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['reaction', 'v_ref_B73'] + L)
    rx = net['rxns'] + ['BIOMASS']
    for a, j in enumerate(nz): w.writerow([rx[j], f'{vref[j]:.6g}'] + [f'{x:.6g}' for x in V['all'][a]])
print(json.dumps({k: summary[k] for k in ('status', 'median_rho', 'beyond_growth', 'mean_abs_rel_dev_from_target_all')}, indent=1))
for k in allm: print(f"{mets[k][:28]:28} none {tab['none'][mets[k]]:+.2f}  all {tab['all'][mets[k]]:+.2f}  held-out {tab['held_out'][mets[k]]:+.2f}")
