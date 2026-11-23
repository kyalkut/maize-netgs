"""Data analyses for the meeting:
 A. G x E signal per BOF metabolite: line consistency WW vs DS + variance partition (line / condition / residual)
 B. Do the model-mapped (BOF) metabolites carry growth information? Ridge CV prediction of DS dry mass."""
import os, sys, json, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
Z = bp.load_zhang(bp.P['zhang'], bp.P['cache']); lines = Z['lines']; L = len(lines)
summ = json.load(open(os.path.join(bp.P['out'], 'run_summary.json')))
mapping = {m['model_metabolite']: m['zhang_id'] for m in bp.read_csv(bp.P['mapping'])}
comp = list(csv.DictReader(open(os.path.join(bp.P['out'], 'bof_components.csv'))))
feat = {c['model_metabolite']: c['zhang_feature'] for c in comp if c['varies'] == 'yes'}
def lv(fid, c):
    v = Z['levels'].get((fid, c)); return None if v is None else v.copy()
# ---------- A. G x E --------------------------------------------------------------------
rowsA = []
for mm, fid in sorted(feat.items()):
    w, d = lv(fid, 'WW'), lv(fid, 'DS'); ok = np.isfinite(w) & np.isfinite(d)
    r = np.corrcoef(w[ok], d[ok])[0, 1]
    X = np.vstack([w[ok], d[ok]]).T                      # lines x 2 conditions, no replicates
    gm = X.mean(); line_eff = X.mean(1) - gm; cond_eff = X.mean(0) - gm
    resid = X - gm - line_eff[:, None] - cond_eff[None, :]
    ss_l = 2 * np.sum(line_eff ** 2); ss_c = X.shape[0] * np.sum(cond_eff ** 2); ss_r = np.sum(resid ** 2); tot = ss_l + ss_c + ss_r
    rowsA.append(dict(metabolite=mm, feature=fid, r_WW_DS=r, pct_line=100 * ss_l / tot, pct_condition=100 * ss_c / tot, pct_residual_GxE=100 * ss_r / tot,
                      drought_shift_log2=float(cond_eff[1] - cond_eff[0])))
# background: all features measured in both conditions
allr = []
for (fid, c), v in Z['levels'].items():
    if c != 'WW': continue
    d = Z['levels'].get((fid, 'DS'))
    if d is None: continue
    ok = np.isfinite(v) & np.isfinite(d)
    if ok.sum() > 300 and np.std(v[ok]) > 0 and np.std(d[ok]) > 0: allr.append(np.corrcoef(v[ok], d[ok])[0, 1])
with open(os.path.join(bp.P['out'], 'analysis_GxE_bof_metabolites.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rowsA[0].keys())); w.writeheader()
    for r in sorted(rowsA, key=lambda r: r['r_WW_DS']): w.writerow({k: (f'{v:.3f}' if isinstance(v, float) else v) for k, v in r.items()})
rA = np.array([r['r_WW_DS'] for r in rowsA])
print(f'A. line consistency WW vs DS: BOF metabolites median r = {np.median(rA):.2f} (range {rA.min():.2f}-{rA.max():.2f}); '
      f'all {len(allr)} features median r = {np.median(allr):.2f}')
for r in sorted(rowsA, key=lambda r: r['r_WW_DS'])[:6]:
    print(f'   low consistency (strong G x E / plasticity): {r["metabolite"]:30} r={r["r_WW_DS"]:.2f}  line {r["pct_line"]:.0f}%  cond {r["pct_condition"]:.0f}%  resid {r["pct_residual_GxE"]:.0f}%')
# ---------- B. prediction of DS dry mass ----------------------------------------------------
dm = np.array([Z['s1'].get(l, {}).get('dry') or np.nan for l in lines], float); ok = np.isfinite(dm)
def matrix(fids, cond='DS'):
    M = []
    for fid in fids:
        v = Z['levels'].get((fid, cond))
        if v is None or np.isfinite(v).mean() < 0.9: continue
        v = v.copy(); v[~np.isfinite(v)] = np.nanmean(v); M.append(v)
    return np.array(M).T
def ridge_cv(X, y, reps=20, folds=5, seed=0):
    rng = np.random.default_rng(seed); X = (X - X.mean(0)) / (X.std(0) + 1e-12); y = y - y.mean(); acc = []
    lam = X.shape[1]  # fixed shrinkage ~ rrBLUP-like (variance ratio ~ number of predictors)
    for _ in range(reps):
        idx = rng.permutation(len(y)); pred = np.zeros(len(y))
        for f in range(folds):
            te = idx[f::folds]; tr = np.setdiff1d(idx, te)
            Xt, yt = X[tr], y[tr]
            beta = np.linalg.solve(Xt.T @ Xt + lam * np.eye(X.shape[1]), Xt.T @ (yt - yt.mean()))
            pred[te] = X[te] @ beta + yt.mean()
        acc.append(np.corrcoef(pred, y)[0, 1])
    return np.mean(acc), np.std(acc)
y = dm[ok]
Xb = matrix(list(feat.values()))[ok]
all_fids = sorted({fid for (fid, c) in Z['levels'] if c == 'DS'}); Xa = matrix(all_fids)[ok]
dr_fids = [f for f in all_fids if Z['dr'].get(f)]; Xd = matrix(dr_fids)[ok]
rng = np.random.default_rng(1); rand = [ridge_cv(Xa[:, rng.choice(Xa.shape[1], Xb.shape[1], replace=False)], y, reps=5, seed=s)[0] for s in range(20)]
res = dict(bof=ridge_cv(Xb, y), all=ridge_cv(Xa, y, reps=5), drought_resp=ridge_cv(Xd, y, reps=5))
print(f'B. ridge 5-fold CV prediction of DS dry mass ({ok.sum()} lines):')
print(f'   {Xb.shape[1]} BOF metabolites:            r = {res["bof"][0]:.3f} +/- {res["bof"][1]:.3f}')
print(f'   {Xb.shape[1]} random metabolites (20x):   r = {np.mean(rand):.3f}  (95% range {np.percentile(rand,2.5):.3f}-{np.percentile(rand,97.5):.3f})')
print(f'   {Xd.shape[1]} drought-responsive:      r = {res["drought_resp"][0]:.3f}')
print(f'   all {Xa.shape[1]} metabolites:           r = {res["all"][0]:.3f}')
json.dump(dict(GxE=dict(median_r_bof=float(np.median(rA)), median_r_all=float(np.median(allr))),
               prediction=dict(bof=res['bof'], random_same_size=[float(np.mean(rand)), float(np.percentile(rand, 2.5)), float(np.percentile(rand, 97.5))],
                               drought_responsive=res['drought_resp'], all=res['all'], n_lines=int(ok.sum()))),
          open(os.path.join(bp.P['out'], 'analysis_summary.json'), 'w'), indent=1)
