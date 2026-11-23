"""Python port of netGS steps 2-4 (Tong, Kueken & Nikoloski 2020; 2_rrBLUP.R, 3_Biomass.m, 4_Correlation.R) on maize DS.

  step 2  rrBLUP per flux of netGS step 1 (all reactions active in B73 + BIOMASS), 50 replicates x 3-fold CV
          (rrBLUP mixed.solve(y, Z = SNPs, K = NULL) == GBLUP with K = ZZ'; REML for the variance ratio, as rrBLUP)
  step 3  per test line: fluxes closest to the PREDICTED fluxes on the B73-active network, with the test line's own BOF
          (as 3_Biomass.m: biofunc = biomreaction(:, test lines)); only reactions predicted at least as well as biomass
          (accuracy in that fold >= the biomass flux's) enter the objective, as the original
            'code' : min | sum_k x_k / pred_k - 1 |          (3_Biomass.m: A row vector -> scalar, as in step 1)
            'paper': min sum_k (x_k / pred_k - 1)^2           (per-reaction relative deviation, as Eq. 2 in step 1)
          s.t. S x = 0, model bounds, biomass >= 0 (+ RuBisCO vc/vo bounds for the '_co' runs). No growth window.
  step 4  netGS accuracy = Pearson r (predicted biomass flux, measured DS dry mass) per fold, mean over 150 folds.
  Classical GS: rrBLUP directly on DS dry mass with the same folds (and on the step-1 biomass flux).

Env: BOF_OUT (variant folder, default output/), NETGS_TAG (step-1 run, default paper_cvxpy_ser), NETGS_OBJ3 (code/paper),
     NETGS_REPS (default 50), NETGS_WORKERS (14). Genotypes: output/snp_kernel_DS.npz (snp_prepare.py).
Outputs in the variant folder: netgs_step2_accuracy_DS_<TAG>.csv (per reaction), netgs_step34_folds_DS_<TAG>_<OBJ3>.csv
(per replicate x fold), netgs_step34_DS_<TAG>_<OBJ3>.json (summary)."""
import os, sys, csv, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
TAG = os.environ.get('NETGS_TAG', 'paper_cvxpy_ser'); OBJ3 = os.environ.get('NETGS_OBJ3', 'code')
R = int(os.environ.get('NETGS_REPS', '50')); NF = 3; OUT = bp.P['out']; MAIN = os.path.join(HERE, 'output')
CO = 'co' in TAG.split('_')
# reactions in the step-3 objective (diagnostics; 'acc' = netGS):
#  acc       step-2 accuracy in this fold >= the biomass flux's (3_Biomass.m; uses the test fold's accuracy)
#  inner     the same rule, accuracy from an inner 3-fold CV on the training lines only (no test information)
#  all       every reaction active in B73;  biomass  the biomass flux only (netGS reduces to GS on the biomass flux)
#  growth    'acc' and |r(flux, biomass flux)| >= 0.95 on the training lines (growth-coupled reactions)
#  nongrowth 'acc' and |r| < 0.95 (reactions with their own variation), plus the biomass flux
SEL = os.environ.get('NETGS_SEL', 'acc'); SUF = '' if SEL == 'acc' else f'_sel-{SEL}'

# ---------------------------------------------------------------- data
def load():
    G = np.load(os.path.join(MAIN, 'snp_kernel_DS.npz'), allow_pickle=True); gl = list(G['lines'])
    rows = list(csv.reader(open(os.path.join(OUT, f'netgs_fluxes_DS_{TAG}.csv'))))
    fl = rows[0][2:]; rx = [r[0] for r in rows[1:]]; vref = np.array([float(r[1]) for r in rows[1:]])
    F = np.array([[float(x) for x in r[2:]] for r in rows[1:]])
    feas = {r['line']: r for r in csv.DictReader(open(os.path.join(OUT, 'feasibility_per_line.csv')))}
    ok = np.all(np.isfinite(F), 0)                           # lines whose step-1 QP failed (no fluxes) are left out
    lines = [l for k, l in enumerate(fl) if l in gl and ok[k]]
    gi = [gl.index(l) for l in lines]; fi = [fl.index(l) for l in lines]
    return dict(K=G['K'][np.ix_(gi, gi)], n_snps=int(G['n_snps']), lines=lines, rx=rx, vref=vref, F=F[:, fi].T,   # lines x reactions
                dm=np.array([float(feas[l]['dry_mass_DS_g']) for l in lines]))

def folds(lines):
    """50 replicates x 3 folds, balanced, fixed seed; saved once and reused by every variant (as netGS foldid.csv)."""
    p = os.path.join(MAIN, 'netgs_foldid_DS.csv')
    if os.path.exists(p):
        rows = list(csv.reader(open(p))); d = {r[0]: [int(x) for x in r[1:]] for r in rows[1:]}
        if all(l in d for l in lines): return np.array([d[l] for l in lines])
    rng = np.random.default_rng(20261120); n = len(lines)
    ids = np.array([rng.permutation(np.arange(n) % NF + 1) for _ in range(50)]).T
    with open(p, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['line'] + [f'rep{r+1}' for r in range(50)])
        for l, row in zip(lines, ids): w.writerow([l] + list(row))
    return ids

# ---------------------------------------------------------------- step 2: rrBLUP (REML, vectorised over traits)
LOGD = np.linspace(-5, 5, 401)
def rrblup(Ktr, Kte, Y):
    """GBLUP with REML (EMMA form, as rrBLUP mixed.solve), one variance ratio per column of Y (n_train x traits)."""
    n = Ktr.shape[0]; Sm = np.eye(n) - 1.0 / n
    xi, U = np.linalg.eigh(Sm @ Ktr @ Sm); xi, U = xi[1:], U[:, 1:]          # n-1 non-zero (X = intercept)
    xi = np.maximum(xi, 0)
    eta2 = (U.T @ Y) ** 2                                                     # (n-1) x T
    d = 10.0 ** LOGD
    ll = np.array([-(n - 1) / 2 * np.log(np.sum(eta2 / (xi[:, None] + di), 0)) - 0.5 * np.sum(np.log(xi + di)) for di in d])
    delta = d[np.argmax(ll, 0)]                                               # per trait
    lam, V = np.linalg.eigh(Ktr); lam = np.maximum(lam, 0)
    a, b = V.T @ Y, V.T @ np.ones(n)
    w = 1.0 / (lam[:, None] + delta[None, :])                                 # n x T
    beta = (b @ (w * a)) / ((b ** 2) @ w)
    hres = V @ (w * (a - b[:, None] * beta[None, :]))                         # H^-1 (y - beta)
    pred = beta[None, :] + Kte @ hres
    const = Y.std(0) < 1e-12 * np.maximum(1, np.abs(Y.mean(0)))
    pred[:, const] = Y.mean(0)[const]
    return pred, delta

def pearson(a, b):
    if np.std(a) < 1e-15 or np.std(b) < 1e-15: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

# ---------------------------------------------------------------- step 3: network QP (one problem per worker, parameters per line)
_W = {}
def set_rx(rx): _W['rx'] = rx
def init_worker():
    import cvxpy as cp, scipy.sparse as sp
    from model_curation import curate
    m = bp.read_sbml(bp.P['model']); curate(m, bp.P['model'], bp.P['gapfill'], bp.P['broken'])
    bp.PARAM['light_uptake'] = 100.0; net = bp.autotrophic_bounds(m)
    rows = list(csv.reader(open(os.path.join(OUT, 'biomreaction_DS.csv')))); bl = rows[0][3:]
    B = np.zeros((m['S'].shape[0], len(bl)))
    for r in rows[1:]:
        B[bp.species_for(m, r[0], (r[1], 'c', 'p', 'm', 'x', 'v'))] -= np.array(r[3:], float)
    for nm, s in (('ATP', -1), ('ADP', 1), ('Orthophosphate', 1)): B[bp.species_for(m, nm, ('c',))] += s * bp.PARAM['atp_maintenance']
    B = B[net['keep']]
    rx = _W['rx']; cols = [net['rxns'].index(r) for r in rx if r != 'BIOMASS']; kb = rx.index('BIOMASS')
    n = len(rx); Sx = np.zeros((net['S'].shape[0], n)); lb = np.zeros(n); ub = np.zeros(n)
    k = [a for a, r in enumerate(rx) if r != 'BIOMASS']
    Sx[:, k] = net['S'][:, cols]; lb[k] = net['lb'][cols]; ub[k] = net['ub'][cols]; lb[kb] = 0; ub[kb] = 1000
    x = cp.Variable(n); bof = cp.Parameter(Sx.shape[0]); w = cp.Parameter(n); msk = cp.Parameter(n)
    cons = [sp.csr_matrix(Sx) @ x + bof * x[kb] == 0, x >= lb, x <= ub]
    if CO:
        jc, jo = rx.index('R_R00024_p'), rx.index('R_R03140_p')
        cons += [x[jc] - 88.0 * x[jo] <= 0, 22.2 * x[jo] - x[jc] <= 0]
    obj = cp.sum_squares(cp.multiply(w, x) - msk) if OBJ3 == 'paper' else cp.abs(w @ x - 1)
    _W.update(cp=cp, prob=cp.Problem(cp.Minimize(obj), cons), x=x, bof=bof, w=w, msk=msk, kb=kb, B=B, bl=bl)

def solve_rep(args):
    r, ids, P, SM, lines = args
    if 'prob' not in _W: init_worker()
    cp = _W['cp']; kb = _W['kb']; out = np.full(len(lines), np.nan); nsel = []; status = []
    for p in range(1, NF + 1):
        te = np.where(ids == p)[0]; sel = SM[p - 1]
        nsel.append(int(sel.sum()))
        for i in te:
            pr = P[i]; ok = sel & (np.abs(pr) > 1e-12)
            wv = np.where(ok, 1.0 / np.where(ok, pr, 1.0), 0.0)
            _W['w'].value = wv; _W['msk'].value = ok.astype(float)
            _W['bof'].value = _W['B'][:, _W['bl'].index(lines[i])]
            try: _W['prob'].solve(solver=cp.CLARABEL)
            except cp.SolverError: pass
            st = _W['prob'].status; status.append(st)
            if st in ('optimal', 'optimal_inaccurate'): out[i] = _W['x'].value[kb]
    return r, out, nsel, status

# ---------------------------------------------------------------- main
if __name__ == '__main__':
    t0 = time.time(); D = load(); lines, F, K, dm, rx = D['lines'], D['F'], D['K'], D['dm'], D['rx']; n = len(lines)
    ids = folds(lines)[:, :R]; kb = rx.index('BIOMASS')
    inref = D['vref'] != 0                                   # serine-pathway extras (v_ref = 0) stay out of the objective, as in step 1
    print(f'{n} lines with genotypes and step-1 fluxes ({D["n_snps"]} SNPs), {len(rx)} reactions, {R} x {NF}-fold CV; TAG {TAG}, step 3 {OBJ3}', flush=True)
    # step 2
    PRED = np.full((R, n, len(rx)), np.nan); ACC = np.full((R, NF, len(rx)), np.nan); SMASK = np.zeros((R, NF, len(rx)), bool)
    PDM = np.full((R, n), np.nan); r_gs = np.full((R, NF), np.nan); r_gsbio = np.full((R, NF), np.nan)
    for r in range(R):
        for p in range(1, NF + 1):
            tr, te = np.where(ids[:, r] != p)[0], np.where(ids[:, r] == p)[0]
            Y = np.column_stack([F[tr], dm[tr]])
            pred, _ = rrblup(K[np.ix_(tr, tr)], K[np.ix_(te, tr)], Y)
            PRED[r][np.ix_(te, np.arange(len(rx)))] = pred[:, :-1]; PDM[r, te] = pred[:, -1]
            ACC[r, p - 1] = [pearson(F[te, k], pred[:, k]) for k in range(len(rx))]
            rb = np.array([abs(pearson(F[tr, k], F[tr, kb])) for k in range(len(rx))])
            acc = ACC[r, p - 1]; base = np.isfinite(acc) & (acc >= acc[kb]) & inref
            if SEL == 'acc': sm = base
            elif SEL == 'all': sm = inref.copy()
            elif SEL == 'biomass': sm = np.zeros(len(rx), bool); sm[kb] = True
            elif SEL == 'growth': sm = base & (rb >= 0.95)
            elif SEL == 'nongrowth': sm = base & ~(rb >= 0.95); sm[kb] = True
            elif SEL == 'inner':
                ia = np.zeros(len(rx)); inner = np.arange(len(tr)) % NF; np.random.default_rng(r * 10 + p).shuffle(inner)
                ip = np.full((len(tr), len(rx)), np.nan)
                for q in range(NF):
                    a_, b_ = tr[inner != q], tr[inner == q]
                    ip[inner == q] = rrblup(K[np.ix_(a_, a_)], K[np.ix_(b_, a_)], F[a_])[0]
                ia = np.array([pearson(F[tr, k], ip[:, k]) for k in range(len(rx))])
                sm = np.isfinite(ia) & (ia >= ia[kb]) & inref
            SMASK[r, p - 1] = sm
            r_gs[r, p - 1] = pearson(dm[te], pred[:, -1]); r_gsbio[r, p - 1] = pearson(dm[te], pred[:, kb])
    print(f'step 2 done {time.time()-t0:.0f}s; classical GS on dry mass r = {np.nanmean(r_gs):.3f}; '
          f'biomass flux predicted r = {np.nanmean(ACC[:, :, kb]):.3f}', flush=True)
    acc_m = np.nanmean(ACC.reshape(-1, len(rx)), 0); acc_s = np.nanstd(ACC.reshape(-1, len(rx)), 0)
    with open(os.path.join(OUT, f'netgs_step2_accuracy_DS_{TAG}.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['reaction', 'v_ref_B73', 'accuracy_mean', 'accuracy_sd', 'share_folds_ge_biomass', 'cv_across_lines', 'r_with_biomass_flux'])
        for k in np.argsort(-np.nan_to_num(acc_m, nan=-9)):
            ge = np.nanmean((ACC[:, :, k] >= ACC[:, :, kb]).astype(float))
            cvl = np.std(F[:, k]) / abs(np.mean(F[:, k])) if np.mean(F[:, k]) != 0 else np.nan
            w.writerow([rx[k], f'{D["vref"][k]:.6g}', f'{acc_m[k]:.4f}', f'{acc_s[k]:.4f}', f'{ge:.3f}', f'{cvl:.4g}', f'{pearson(F[:, k], F[:, kb]):.4f}'])
    # step 3 (parallel over replicates)
    _W['rx'] = rx
    from multiprocessing import Pool
    BIO = np.full((R, n), np.nan); NSEL = np.zeros((R, NF), int); ST = []
    jobs = [(r, ids[:, r], PRED[r], SMASK[r], lines) for r in range(R)]
    with Pool(int(os.environ.get('NETGS_WORKERS', '14')), initializer=set_rx, initargs=(rx,)) as pool:
        for r, out, nsel, status in pool.imap_unordered(solve_rep, jobs):
            BIO[r] = out; NSEL[r] = nsel; ST += status
            print(f'  rep {r+1} done {time.time()-t0:.0f}s', flush=True)
    # step 4
    r_net = np.full((R, NF), np.nan)
    for r in range(R):
        for p in range(1, NF + 1):
            te = np.where(ids[:, r] == p)[0]; ok = np.isfinite(BIO[r, te])
            if ok.sum() > 3: r_net[r, p - 1] = pearson(dm[te][ok], BIO[r, te][ok])
    with open(os.path.join(OUT, f'netgs_step34_folds_DS_{TAG}_{OBJ3}{SUF}.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['replicate', 'fold', 'r_netgs', 'r_classical_gs_drymass', 'r_gs_biomass_flux', 'acc_biomass_flux', 'n_reactions_selected'])
        for r in range(R):
            for p in range(NF):
                w.writerow([r + 1, p + 1, f'{r_net[r, p]:.4f}', f'{r_gs[r, p]:.4f}', f'{r_gsbio[r, p]:.4f}', f'{ACC[r, p, kb]:.4f}', NSEL[r, p]])
    S = dict(tag=TAG, step3_objective=OBJ3, selection=SEL, n_lines=n, n_snps=D['n_snps'], n_reactions=len(rx), replicates=R, folds=NF,
             netgs_accuracy_mean=float(np.nanmean(r_net)), netgs_accuracy_sd=float(np.nanstd(r_net)),
             classical_gs_drymass_mean=float(np.nanmean(r_gs)), classical_gs_drymass_sd=float(np.nanstd(r_gs)),
             gs_on_biomass_flux_mean=float(np.nanmean(r_gsbio)), biomass_flux_prediction_accuracy_mean=float(np.nanmean(ACC[:, :, kb])),
             netgs_minus_gs_mean=float(np.nanmean(r_net - r_gs)), share_folds_netgs_better=float(np.nanmean((r_net > r_gs).astype(float))),
             n_reactions_selected_median=float(np.median(NSEL)), n_reactions_selected_range=[int(NSEL.min()), int(NSEL.max())],
             step2_accuracy_median_over_reactions=float(np.nanmedian(acc_m)),
             qp_status={s: ST.count(s) for s in set(ST)}, seconds=round(time.time() - t0))
    json.dump(S, open(os.path.join(OUT, f'netgs_step34_DS_{TAG}_{OBJ3}{SUF}.json'), 'w'), indent=1)
    print(json.dumps(S, indent=1))
