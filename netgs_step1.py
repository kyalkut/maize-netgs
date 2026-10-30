"""
Python port of netGS step 1 (Tong, Kueken & Nikoloski 2020; 1_FluxDist.m): genotype-specific flux distributions.

  1a reference flux v_ref: FBA (max biomass) for the reference genotype; keep non-zero reactions
  1b s_model: per genotype, max biomass with its own BOF  -> mean over genotypes (netGS) or min
  1c growth rescaling: target_i = (s_model - delta) * M_i / max(M);  window +/- 1.645*SEM (90%)
  1d per genotype: stay close to v_ref subject to S_i v = 0, bounds, biomass window, extra constraints

Two objectives for 1d:
  'paper' : min sum_k ((v_k - v_ref,k) / v_ref,k)^2      (relative squared deviation per reaction; Methods Eq. 2)
  'code'  : min | sum_k v_k / v_ref,k  - 1 |               (what 1_FluxDist.m computes: A = 1./fluxc.' is a 1xn ROW,
                                                             so A*x is a scalar and norm(A*x - b) is |sum - 1|)
"""
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
from qp_admm import solve_qp

def fba(S, lb, ub, c, A_ub=None, b_ub=None, pfba=True):
    res = linprog(-c, A_eq=S, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub, bounds=list(zip(lb, ub)), method='highs')
    if res.status != 0: return None, 0.0
    zmax = -res.fun
    if not pfba: return res.x, zmax
    # parsimonious FBA: fix biomass at optimum, minimise total absolute flux (unique, like 1_FluxDist_pFBA.m)
    n = S.shape[1]; jb = int(np.argmax(c))
    A = sp.hstack([sp.csr_matrix(S), -sp.csr_matrix(S)]).tocsr()
    lbp = np.maximum(lb, 0); ubp = np.maximum(ub, 0); lbn = np.maximum(-ub, 0); ubn = np.maximum(-lb, 0)
    lbp[jb] = ubp[jb] = zmax * 0.999999; lbn[jb] = ubn[jb] = 0
    Aub = None if A_ub is None else sp.hstack([sp.csr_matrix(A_ub), -sp.csr_matrix(A_ub)]).tocsr()
    r2 = linprog(np.ones(2 * n), A_eq=A, b_eq=np.zeros(S.shape[0]), A_ub=Aub, b_ub=b_ub,
                 bounds=list(zip(lbp, ubp)) + list(zip(lbn, ubn)), method='highs')
    return (r2.x[:n] - r2.x[n:], zmax) if r2.status == 0 else (res.x, zmax)

def genotype_fluxes(S, lb, ub, jbio, bofs, measured, ref, extra_ub=None, objective='paper', delta=1.1e-4,
                    ref_bof=None, smodel='mean', verbose=False, qp_subset=None):
    """S: metabolites x reactions (biomass column jbio is replaced per genotype by bofs[:, i]).
    extra_ub: (A, b) additional inequalities A v <= b (e.g. C3 ratio bounds). Returns dict."""
    n = S.shape[1]; c = np.zeros(n); c[jbio] = 1
    A_ub, b_ub = extra_ub if extra_ub is not None else (None, None)
    S0 = S.copy();
    if ref_bof is not None: S0[:, jbio] = ref_bof
    vref, _ = fba(S0, lb, ub, c, A_ub, b_ub)
    nz = np.where(np.abs(vref) > 1e-5)[0]
    # 1b: model maximum per genotype
    zmax = []
    for i in range(bofs.shape[1]):
        Si = S.copy(); Si[:, jbio] = bofs[:, i]
        r = linprog(-c, A_eq=Si, b_eq=np.zeros(S.shape[0]), A_ub=A_ub, b_ub=b_ub, bounds=list(zip(lb, ub)), method='highs')
        zmax.append(-r.fun if r.status == 0 else 0.0)
    zmax = np.array(zmax); s_model = zmax.mean() if smodel == 'mean' else zmax.min()
    # 1c: growth rescaling (two values: s_model and s_measurement = max measured)
    target = (s_model - delta) * measured / np.nanmax(measured)
    sem = np.nanstd(target, ddof=1) / np.sqrt(np.sum(np.isfinite(target))); half = 1.645 * sem
    # 1d: per genotype, restricted to reactions active in the reference
    Sn = S[:, nz]; lbn, ubn = lb[nz], ub[nz]; vr = vref[nz]; kb = int(np.where(nz == jbio)[0][0])
    Aun = None if A_ub is None else A_ub[:, nz]
    out = np.full((len(nz), bofs.shape[1]), np.nan); info = []
    todo = range(bofs.shape[1]) if qp_subset is None else qp_subset
    for i in range(bofs.shape[1]):
        if i not in todo: info.append('skipped'); continue
        if not np.isfinite(target[i]): info.append('no growth data'); continue
        Si = Sn.copy(); Si[:, kb] = bofs[:, i]
        lo_b, hi_b = target[i] - half, target[i] + half
        lbi, ubi = lbn.copy(), ubn.copy(); lbi[kb] = max(lbi[kb], lo_b); ubi[kb] = min(ubi[kb], hi_b)
        if objective == 'paper':
            # well-conditioned form: y = v / v_ref (relative flux); min 1/2 ||y - 1||^2
            D = sp.diags(vr)
            rows = [sp.csr_matrix(Si) @ D, sp.eye(len(nz))]
            ylo = np.where(vr > 0, lbi / vr, ubi / vr); yhi = np.where(vr > 0, ubi / vr, lbi / vr)
            l = [np.zeros(Si.shape[0]), ylo]; u = [np.zeros(Si.shape[0]), yhi]
            if Aun is not None: rows.append(sp.csr_matrix(Aun) @ D); l.append(np.full(len(b_ub), -np.inf)); u.append(b_ub)
            y, st = solve_qp(sp.eye(len(nz)), -np.ones(len(nz)), sp.vstack(rows), np.concatenate(l), np.concatenate(u),
                             x0=np.ones(len(nz)), rho=1.0, max_iter=40000, eps_abs=1e-7, eps_rel=1e-6)
            x = y * vr
            info.append(st['status'])
        else:  # 'code': min |a.x - 1| as LP with auxiliary t
            a = 1.0 / vr; m = len(nz)
            cc = np.r_[np.zeros(m), 1.0]
            Aeq = np.hstack([Si, np.zeros((Si.shape[0], 1))])
            Aub2 = np.vstack([np.r_[a, -1.0], np.r_[-a, -1.0]]); bub2 = np.array([1.0, -1.0])
            if Aun is not None:
                Aub2 = np.vstack([Aub2, np.hstack([Aun, np.zeros((Aun.shape[0], 1))])]); bub2 = np.r_[bub2, b_ub]
            r = linprog(cc, A_eq=Aeq, b_eq=np.zeros(Si.shape[0]), A_ub=Aub2, b_ub=bub2,
                        bounds=list(zip(lbi, ubi)) + [(0, None)], method='highs')
            x = r.x[:m] if r.status == 0 else np.full(m, np.nan); info.append('solved' if r.status == 0 else 'infeasible')
        out[:, i] = x
    return dict(vref=vref, nonzero=nz, zmax=zmax, s_model=s_model, target=target, halfwidth=half, fluxes=out, status=info)
