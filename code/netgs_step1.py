"""
Python port of netGS step 1 (Tong, Kueken & Nikoloski 2020; 1_FluxDist.m): genotype-specific flux distributions.

  1a reference flux v_ref: FBA (max biomass) for the reference genotype; keep non-zero reactions
  1b s_model: per genotype, max biomass with its own BOF  -> mean over genotypes (netGS) or min
  1c growth rescaling: target_i = (s_model - delta) * M_i / max(M);  window +/- 1.645*SEM (90%)
  1d per genotype: stay close to v_ref subject to S_i v = 0, bounds, biomass window, extra constraints

QP solver for 'paper': solver='osqp' (default; pip install osqp) or 'admm' (own solver in qp_admm.py).
Solver for 'code': solver='cvxpy' (as the original: CVX-style norm, interior point; pip install cvxpy) or HiGHS LP (any other value).

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

def solve_qp_osqp(P, q, C, l, u, x0=None):
    """Same QP as qp_admm.solve_qp (min 1/2 x'Px + q'x, l <= Cx <= u), solved with OSQP + solution polishing."""
    import osqp
    prob = osqp.OSQP()
    opts = dict(eps_abs=1e-8, eps_rel=1e-8, max_iter=200000, verbose=False)
    try: prob.setup(sp.csc_matrix(P), q, sp.csc_matrix(C), l, u, polishing=True, **opts)   # osqp >= 1.0
    except TypeError: prob.setup(sp.csc_matrix(P), q, sp.csc_matrix(C), l, u, polish=True, **opts)
    if x0 is not None: prob.warm_start(x=x0)
    res = prob.solve()
    st = res.info.status
    return res.x, dict(status='solved' if st == 'solved' else st, iter=res.info.iter,
                       r_prim=res.info.prim_res, r_dual=res.info.dual_res)

def genotype_fluxes(S, lb, ub, jbio, bofs, measured, ref, extra_ub=None, objective='paper', delta=1.1e-4,
                    ref_bof=None, smodel='mean', verbose=False, qp_subset=None, solver='osqp', extra_rxns=()):
    """S: metabolites x reactions (biomass column jbio is replaced per genotype by bofs[:, i]).
    extra_ub: (A, b) additional inequalities A v <= b (e.g. C3 ratio bounds).
    extra_rxns: reaction indices added to the reference-active set although v_ref = 0 (e.g. the phosphorylated serine
      pathway, which high-purine lines need). 'paper': penalised as (v / c)^2, c = median |v_ref| of active reactions;
      'code': left out of the sum (A = 0), since 1/v_ref is undefined. Returns dict."""
    n = S.shape[1]; c = np.zeros(n); c[jbio] = 1
    A_ub, b_ub = extra_ub if extra_ub is not None else (None, None)
    S0 = S.copy();
    if ref_bof is not None: S0[:, jbio] = ref_bof
    vref, _ = fba(S0, lb, ub, c, A_ub, b_ub)
    nz = np.where(np.abs(vref) > 1e-5)[0]
    cscale = np.median(np.abs(vref[nz]))
    ext = np.array(sorted(set(extra_rxns) - set(nz)), dtype=int); nz = np.sort(np.r_[nz, ext])
    isx = np.isin(nz, ext)
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
            # well-conditioned form: y = v / v_ref (relative flux); min 1/2 ||y - 1||^2  (extra reactions: y = v / c, target 0)
            sc = np.where(isx, cscale, vr); q = np.where(isx, 0.0, -1.0)
            D = sp.diags(sc)
            rows = [sp.csr_matrix(Si) @ D, sp.eye(len(nz))]
            ylo = np.where(sc > 0, lbi / sc, ubi / sc); yhi = np.where(sc > 0, ubi / sc, lbi / sc)
            l = [np.zeros(Si.shape[0]), ylo]; u = [np.zeros(Si.shape[0]), yhi]
            if Aun is not None: rows.append(sp.csr_matrix(Aun) @ D); l.append(np.full(len(b_ub), -np.inf)); u.append(b_ub)
            if solver == 'cvxpy':   # interior point (Clarabel), same QP
                import cvxpy as cp
                y_ = cp.Variable(len(nz)); Cm = sp.vstack(rows).tocsr(); lo_, hi_ = np.concatenate(l), np.concatenate(u)
                fin_lo, fin_hi = np.isfinite(lo_), np.isfinite(hi_)
                pr = cp.Problem(cp.Minimize(0.5 * cp.sum_squares(y_ + q)), [Cm[fin_lo] @ y_ >= lo_[fin_lo], Cm[fin_hi] @ y_ <= hi_[fin_hi]])
                try: pr.solve(solver=cp.CLARABEL)
                except cp.SolverError: pass
                ok = pr.status in ('optimal', 'optimal_inaccurate')
                y = np.asarray(y_.value) if ok else np.full(len(nz), np.nan)
                st = dict(status='solved' if pr.status == 'optimal' else pr.status)
            elif solver == 'osqp':
                y, st = solve_qp_osqp(sp.eye(len(nz)), q, sp.vstack(rows), np.concatenate(l), np.concatenate(u),
                                      x0=-q)
            else:
                y, st = solve_qp(sp.eye(len(nz)), q, sp.vstack(rows), np.concatenate(l), np.concatenate(u),
                                 x0=-q, rho=1.0, max_iter=40000, eps_abs=1e-7, eps_rel=1e-6)
            x = y * sc
            info.append(st['status'])
        elif solver == 'cvxpy':  # 'code' exactly as 1_FluxDist.m: CVX minimize(norm(A*x-b)), interior point (Clarabel ~ SDPT3)
            import cvxpy as cp
            a = np.where(isx, 0.0, 1.0 / np.where(isx, 1.0, vr)); x_ = cp.Variable(len(nz))
            cons = [sp.csr_matrix(Si) @ x_ == 0, x_ >= lbi, x_ <= ubi]
            if Aun is not None: cons.append(Aun @ x_ <= b_ub)
            pr = cp.Problem(cp.Minimize(cp.norm(a @ x_ - 1)), cons)
            try: pr.solve(solver=cp.CLARABEL)
            except cp.SolverError: pass
            ok = pr.status in ('optimal', 'optimal_inaccurate')
            x = np.asarray(x_.value) if ok else np.full(len(nz), np.nan); info.append('solved' if pr.status == 'optimal' else pr.status)
        else:  # 'code': min |a.x - 1| as LP with auxiliary t (HiGHS simplex -> a vertex of the optimal set)
            a = np.where(isx, 0.0, 1.0 / np.where(isx, 1.0, vr)); m = len(nz)
            cc = np.r_[np.zeros(m), 1.0]
            Aeq = np.hstack([Si, np.zeros((Si.shape[0], 1))])
            Aub2 = np.vstack([np.r_[a, -1.0], np.r_[-a, -1.0]]); bub2 = np.array([1.0, -1.0])
            if Aun is not None:
                Aub2 = np.vstack([Aub2, np.hstack([Aun, np.zeros((Aun.shape[0], 1))])]); bub2 = np.r_[bub2, b_ub]
            r = linprog(cc, A_eq=Aeq, b_eq=np.zeros(Si.shape[0]), A_ub=Aub2, b_ub=bub2,
                        bounds=list(zip(lbi, ubi)) + [(0, None)], method='highs')
            x = r.x[:m] if r.status == 0 else np.full(m, np.nan); info.append('solved' if r.status == 0 else 'infeasible')
        out[:, i] = x
    return dict(vref=vref, nonzero=nz, extra=ext, zmax=zmax, s_model=s_model, target=target, halfwidth=half, fluxes=out, status=info)
