"""Merge chunked netGS step-1 runs, write the flux table, summarise and plot.
NETGS_TAG picks the run (default paper_cvxpy_ser; part files from run_netgs_maize.py with NETGS_CHUNK)."""
import os, glob, json, csv, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = os.environ.get('BOF_OUT') or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
TAG = os.environ.get('NETGS_TAG', 'paper_cvxpy_ser')
parts = sorted(glob.glob(os.path.join(OUT, f'netgs_part_DS_{TAG}_[0-9]*.npy')))
P0 = np.load(parts[0], allow_pickle=True).item(); F = np.full_like(P0['fluxes'], np.nan); status = list(P0['status'])
for p in parts:
    d = np.load(p, allow_pickle=True).item()
    for i, s in enumerate(d['status']):
        if s != 'skipped': F[:, i] = d['fluxes'][:, i]; status[i] = s
nz, vref, keep, lines, rx = P0['nonzero'], P0['vref'], P0['keep'], P0['lines'], P0['rx']
names = [rx[j] for j in nz]; kb = names.index('BIOMASS'); L = [lines[k] for k in keep]
with open(os.path.join(OUT, f'netgs_fluxes_DS_{TAG}.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['reaction', 'v_ref_B73'] + L)
    for a, j in enumerate(nz): w.writerow([rx[j], f'{vref[j]:.6g}'] + [f'{x:.6g}' for x in F[a]])
bio = F[kb]; tgt = P0['target']; half = P0['half']
ok_window = np.mean(np.abs(bio - tgt) <= half * 1.001)
# variation beyond growth: normalise each line's fluxes by its biomass flux
Fn = F / bio[None, :]
cv_raw = np.nanstd(F, 1) / np.abs(np.nanmean(F, 1)); cv_norm = np.nanstd(Fn, 1) / np.abs(np.nanmean(Fn, 1))
r_bio = np.array([np.corrcoef(F[a], bio)[0, 1] if np.nanstd(F[a]) > 0 else np.nan for a in range(len(nz))])
order = np.argsort(-np.nan_to_num(cv_norm))
summary = dict(n_lines=int(len(L)), n_reactions=int(len(nz)), status={s: status.count(s) for s in set(status)},
               biomass_within_window=float(ok_window), s_model=float(P0['s_model']), halfwidth=float(half),
               median_cv_raw=float(np.nanmedian(cv_raw)), median_cv_after_growth_normalisation=float(np.nanmedian(cv_norm)),
               share_reactions_r_bio_gt_0_9=float(np.nanmean(r_bio > 0.9)),
               top_variable_reactions_after_normalisation=[(names[a], float(cv_norm[a])) for a in order[:12]])
json.dump(summary, open(os.path.join(OUT, f'netgs_summary_DS_{TAG}.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != 'top_variable_reactions_after_normalisation'}, indent=1))
for n, v in summary['top_variable_reactions_after_normalisation']: print(f'  {n:22} CV {v:.3f}')
# figure
BLUE, ORANGE, INK, INK2, GRID, SURF = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF, 'axes.facecolor': SURF,
                     'axes.edgecolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2})
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
ax = axs[0]; dm = np.array([float(x) for x in tgt]) / (P0['s_model'] - 1.1e-4)
ax.scatter(tgt, bio, s=10, color=BLUE, alpha=0.7); lo, hi = np.nanmin(tgt), np.nanmax(tgt)
ax.plot([lo, hi], [lo, hi], color=INK2, lw=1, ls='--'); ax.set_xlabel('target biomass flux (rescaled dry mass)'); ax.set_ylabel('biomass flux in solution')
ax.set_title('Every line meets its growth target', fontsize=10, loc='left', color=INK); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
ax = axs[1]; v = np.sort(cv_norm[np.isfinite(cv_norm)])[::-1]
ax.plot(np.arange(1, len(v) + 1), v, color=ORANGE, lw=2)
ax.set_xlabel('reactions, ranked'); ax.set_ylabel('CV across lines (fluxes per unit growth)')
ax.set_title('Flux variation beyond growth differences', fontsize=10, loc='left', color=INK); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
fig.tight_layout(); FIG = 'fig_netgs_fluxes.png' if TAG == 'paper_cvxpy_ser' else f'fig_netgs_fluxes_{TAG}.png'   # main run keeps the plain name
fig.savefig(os.path.join(OUT, FIG), dpi=160)
print(f'written netgs_fluxes_DS_{TAG}.csv, netgs_summary_DS_{TAG}.json, {FIG}')
