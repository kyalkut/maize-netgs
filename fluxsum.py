"""Flux-sums of the measured metabolites from netGS step 1 fluxes (Seyis, Razaghi-Moghadam & Nikoloski 2025).
  phi_i = 1/2 * sum_j |S_ij * v_j|, summed over all compartments of metabolite i (turnover through its pools)
Per line and DS run, then across lines: phi vs the measured DS log2 level (Spearman), raw and per unit growth.
Main variant excludes the biomass reaction (its coefficients come from the same data -> circular); 'with_bof' includes it.
NETGS_TAG picks the flux file (default paper_cvxpy_ser)."""
import os, sys, csv, json, numpy as np
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
from model_curation import curate
TAG = os.environ.get('NETGS_TAG', 'paper_cvxpy_ser'); OUT = bp.P['out']
m = bp.read_sbml(bp.P['model']); curate(m, bp.P['model'], bp.P['gapfill'], bp.P['broken'])
bp.PARAM['light_uptake'] = 100.0; net = bp.autotrophic_bounds(m)
base = np.array(m['base'])[net['keep']]
# fluxes: reactions x lines
rows = list(csv.reader(open(os.path.join(OUT, f'netgs_fluxes_DS_{TAG}.csv'))))
lines = rows[0][2:]; rx = [r[0] for r in rows[1:]]; F = np.array([[float(x) for x in r[2:]] for r in rows[1:]])
kb = rx.index('BIOMASS'); bio = F[kb]
cols = [net['rxns'].index(r) for r in rx if r != 'BIOMASS']; Fm = F[[a for a, r in enumerate(rx) if r != 'BIOMASS']]
Sn = net['S'][:, cols]
# biomass column per line (consumption of BOF components; ATP maintenance left out of phi)
bofrows = list(csv.reader(open(os.path.join(OUT, 'biomreaction_DS.csv')))); bl = bofrows[0][3:]
bofcoef = {r[0]: np.array([float(x) for x in r[3:]])[[bl.index(l) for l in lines]] for r in bofrows[1:]}
# measured metabolites: all mapped Zhang features; mark the 31 that enter the BOF
Z = bp.load_zhang(bp.P['zhang'], bp.P['cache']); zl = Z['lines']; li = [zl.index(l) for l in lines]
comp = {c['model_metabolite']: c for c in csv.DictReader(open(os.path.join(OUT, 'bof_components.csv')))}
res = []
for mp in bp.read_csv(bp.P['mapping']):
    name, fid = mp['model_metabolite'], mp['zhang_id']
    rws = np.where(base == name)[0]
    lev = Z['levels'].get((fid, 'DS'))
    if len(rws) == 0 or lev is None: continue
    lev = lev[li]
    phi = 0.5 * np.abs(Sn[rws][:, :, None] * Fm[None, :, :]).sum(axis=(0, 1))       # lines
    phi_bof = phi + 0.5 * bofcoef.get(name, np.zeros(len(lines))) * bio
    ok = np.isfinite(lev)
    def rho(x):
        if np.nanstd(x[ok]) < 1e-12 * max(1e-30, np.nanmean(np.abs(x[ok]))): return np.nan
        return spearmanr(x[ok], lev[ok])[0]
    in_bof = comp.get(name, {}).get('varies') == 'yes' and comp[name].get('zhang_feature') == fid
    res.append(dict(metabolite=name, feature=fid, in_bof='yes' if in_bof else 'no', n_reactions=int((np.abs(Sn[rws]) > 0).any(0).sum()),
                    phi_median=float(np.median(phi)), phi_cv=float(np.std(phi) / np.mean(phi)) if np.mean(phi) > 0 else np.nan,
                    rho_level_phi=rho(phi), rho_level_phi_per_growth=rho(phi / bio), rho_level_phi_with_bof=rho(phi_bof),
                    rho_level_growth=spearmanr(bio[ok], lev[ok])[0], n_lines=int(ok.sum()), _phi=phi))
# per-line table of phi (metabolite x line)
with open(os.path.join(OUT, f'fluxsum_DS_{TAG}.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['metabolite', 'zhang_feature', 'in_bof'] + lines)
    for r in res: w.writerow([r['metabolite'], r['feature'], r['in_bof']] + [f'{x:.6g}' for x in r['_phi']])
keys = [k for k in res[0] if not k.startswith('_')]
with open(os.path.join(OUT, f'fluxsum_vs_level_DS_{TAG}.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in sorted(res, key=lambda r: -np.nan_to_num(r['rho_level_phi'], nan=-9)):
        w.writerow({k: (f'{r[k]:.3f}' if isinstance(r[k], float) else r[k]) for k in keys})
def summ(sel):
    out = {}
    for k in ('rho_level_phi', 'rho_level_phi_per_growth', 'rho_level_phi_with_bof', 'rho_level_growth'):
        v = np.array([r[k] for r in sel], float); v = v[np.isfinite(v)]
        out[k] = dict(n=int(len(v)), median=float(np.median(v)) if len(v) else None, n_pos_gt_0_2=int(np.sum(v > 0.2)), n_neg_lt_m0_2=int(np.sum(v < -0.2)))
    return out
S_ = dict(tag=TAG, n_lines=len(lines), n_metabolites=len(res), n_constant_phi=int(sum(1 for r in res if not np.isfinite(r['rho_level_phi']))),
          bof=summ([r for r in res if r['in_bof'] == 'yes']), not_bof=summ([r for r in res if r['in_bof'] == 'no']))
json.dump(S_, open(os.path.join(OUT, f'fluxsum_summary_DS_{TAG}.json'), 'w'), indent=1)
print(json.dumps(S_, indent=1))
for r in sorted(res, key=lambda r: -np.nan_to_num(r['rho_level_phi'], nan=-9)):
    print(f"{r['metabolite'][:28]:28} bof={r['in_bof']:3} rxns={r['n_reactions']:2} phi={r['phi_median']:.3g} cv={r['phi_cv']:.3f} "
          f"rho={r['rho_level_phi']:.2f} per_growth={r['rho_level_phi_per_growth']:.2f} with_bof={r['rho_level_phi_with_bof']:.2f} growth={r['rho_level_growth']:.2f}")
# figure: per metabolite, rho(level, phi) vs number of reactions touching it (low-turnover metabolites echo their BOF coefficient)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BLUE, ORANGE, INK, INK2, GRID, SURF = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF, 'axes.facecolor': SURF,
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
for ax, key, title in ((axs[0], 'rho_level_phi', 'Flux-sum vs measured level'), (axs[1], 'rho_level_phi_per_growth', 'Flux-sum per unit growth vs measured level')):
    for sel, col, lab in (('yes', BLUE, 'in BOF (31)'), ('no', ORANGE, 'measured, not in BOF')):
        pts = [(r['n_reactions'], r[key]) for r in res if r['in_bof'] == sel and np.isfinite(r[key])]
        if pts: ax.scatter(*zip(*pts), s=18, color=col, alpha=0.8, label=lab)
    ax.set_xscale('symlog', linthresh=1); ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel('reactions carrying the metabolite (active network)'); ax.set_ylabel('Spearman r across 371 DS lines')
    ax.set_title(title, fontsize=10, loc='left', color=INK); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_ylim(-1.05, 1.05)
axs[0].legend(frameon=False, fontsize=8, loc='lower right')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig_fluxsum.png'), dpi=160)
