"""Small OSQP-style ADMM solver for convex QPs:  min 1/2 x'Px + q'x  s.t.  l <= Cx <= u
with adaptive rho (as in OSQP). Only numpy/scipy are needed. Used for the netGS per-line quadratic programs."""
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

def _factor(P, C, sigma, rho_v):
    n = P.shape[0]
    K = sp.bmat([[P + sigma * sp.eye(n), C.T], [C, -sp.diags(1.0 / rho_v)]], format='csc')
    return splu(K)

def solve_qp(P, q, C, l, u, rho=0.1, sigma=1e-6, alpha=1.6, max_iter=20000, eps_abs=1e-7, eps_rel=1e-7,
             x0=None, adapt_every=100):
    P = sp.csc_matrix(P); C = sp.csc_matrix(C); n = P.shape[0]; m = C.shape[0]
    eq = np.isclose(l, u)
    rho_v = np.where(eq, 1e3 * rho, rho); lu = _factor(P, C, sigma, rho_v)
    x = np.zeros(n) if x0 is None else x0.copy(); z = np.clip(C @ x, l, u); y = np.zeros(m)
    for k in range(1, max_iter + 1):
        rhs = np.concatenate([sigma * x - q, z - y / rho_v])
        sol = lu.solve(rhs); xt = sol[:n]; nu = sol[n:]
        zt = z + (nu - y) / rho_v
        x_new = alpha * xt + (1 - alpha) * x
        z_new = np.clip(alpha * zt + (1 - alpha) * z + y / rho_v, l, u)
        y = y + rho_v * (alpha * zt + (1 - alpha) * z - z_new)
        x, z = x_new, z_new
        if k % 25 == 0 or k == max_iter:
            Cx = C @ x; Px = P @ x; Cty = C.T @ y
            r_p = np.max(np.abs(Cx - z)); r_d = np.max(np.abs(Px + q + Cty))
            sc_p = max(np.max(np.abs(Cx)), np.max(np.abs(z)), 1e-12)
            sc_d = max(np.max(np.abs(Px)), np.max(np.abs(Cty)), np.max(np.abs(q)), 1e-12)
            if r_p <= eps_abs + eps_rel * sc_p and r_d <= eps_abs + eps_rel * sc_d:
                return x, dict(status='solved', iter=k, r_prim=r_p, r_dual=r_d)
            if k % adapt_every == 0:
                ratio = np.sqrt((r_p / sc_p) / max(r_d / sc_d, 1e-16))
                if ratio > 5 or ratio < 0.2:
                    rho_v = np.clip(rho_v * ratio, 1e-6, 1e6); lu = _factor(P, C, sigma, rho_v)
    return x, dict(status='max_iter', iter=max_iter, r_prim=r_p, r_dual=r_d)
