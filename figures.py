"""Meeting figures (static PNG). Palette: reference categorical slots (blue #2a78d6, orange #eb6834), light surface."""
import os, json, csv, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
BLUE, ORANGE, INK, INK2, GRID, SURF = '#2a78d6', '#eb6834', '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'figure.facecolor': SURF, 'axes.facecolor': SURF})
# --- Fig A: G x E consistency
rows = list(csv.DictReader(open(os.path.join(OUT, 'analysis_GxE_bof_metabolites.csv'))))
s = json.load(open(os.path.join(OUT, 'analysis_summary.json')))
rows.sort(key=lambda r: float(r['r_WW_DS']))
fig, ax = plt.subplots(figsize=(7.6, 7.2))
y = np.arange(len(rows)); r = np.array([float(x['r_WW_DS']) for x in rows])
aa = np.array([x['metabolite'].startswith('L-') for x in rows])
ax.hlines(y, 0, r, color=GRID, lw=2)
ax.scatter(r[~aa], y[~aa], s=40, color=BLUE, zorder=3, label='other BOF metabolites')
ax.scatter(r[aa], y[aa], s=40, color=ORANGE, zorder=3, label='amino acids (Table III)')
ax.axvline(s['GxE']['median_r_all'], color=INK2, ls='--', lw=1)
ax.text(s['GxE']['median_r_all'] + 0.01, 3.5, f"median of all\n3,116 metabolites\n(r = {s['GxE']['median_r_all']:.2f})", color=INK2, fontsize=8, va='center')
ax.set_yticks(y); ax.set_yticklabels([x['metabolite'] for x in rows], fontsize=8)
ax.set_xlabel('correlation of line values, WW vs DS\n(low = lines change their ranking under drought)')
ax.set_xlim(0, 1); ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.legend(loc='upper right', bbox_to_anchor=(1.0, 0.62), frameon=False, fontsize=8)
ax.set_title('Line consistency WW vs DS, 31 BOF metabolites', fontsize=10, color=INK, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig_GxE_consistency.png'), dpi=160); plt.close(fig)
# --- Fig B: prediction of DS dry mass
p = s['prediction']
labels = ['31 BOF metabolites', '31 random metabolites', '825 drought-responsive', 'all 3,134 metabolites']
vals = [p['bof'][0], p['random_same_size'][0], p['drought_responsive'][0], p['all'][0]]
fig, ax = plt.subplots(figsize=(6.4, 2.8))
yy = np.arange(len(labels))[::-1]
cols = [ORANGE, BLUE, BLUE, BLUE]
ax.barh(yy, vals, color=cols, height=0.55)
ax.errorbar(p['random_same_size'][0], yy[1], xerr=[[p['random_same_size'][0] - p['random_same_size'][1]], [p['random_same_size'][2] - p['random_same_size'][0]]],
            color=INK, capsize=3, lw=1)
ends = [vals[0], p['random_same_size'][2], vals[2], vals[3]]
for v, e, yv in zip(vals, ends, yy): ax.text(e + 0.012, yv, f'{v:.2f}', va='center', color=INK, fontsize=8)
ax.set_yticks(yy); ax.set_yticklabels(labels); ax.set_xlim(0, 0.65)
ax.set_xlabel('prediction accuracy of DS dry mass (r, ridge, 5-fold CV)')
ax.grid(axis='x', color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title('Growth information in the metabolites the model can use', fontsize=10, color=INK, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig_prediction_drymass.png'), dpi=160); plt.close(fig)
print('figures written')
