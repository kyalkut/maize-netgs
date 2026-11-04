"""netGS step 1 on maize (DS): curated C4GEM (single cell), line-specific DS BOFs, DS dry mass, reference B73."""
import sys, os, csv, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
from model_curation import curate
from netgs_step1 import genotype_fluxes
OBJ = os.environ.get('NETGS_OBJ', 'paper'); NMAX = int(os.environ.get('NETGS_N', '100000'))
cond = os.environ.get('NETGS_COND', 'DS')
SOLVER = os.environ.get('NETGS_SOLVER', 'cvxpy')
TAG = OBJ + ('_osqp' if SOLVER == 'osqp' and OBJ == 'paper' else '_cvxpy' if SOLVER == 'cvxpy' else '')
# phosphorylated serine pathway: zero in B73 but needed by 67 high-purine lines on the B73-active network ('none' = off)
EXTRA = os.environ.get('NETGS_EXTRA', 'R_R01513_c,R_R04173_c,R_R00582_c')
EXTRA = [] if EXTRA.lower() == 'none' else EXTRA.split(',')
if EXTRA: TAG += '_ser'
m = bp.read_sbml(bp.P['model']); curate(m, bp.P['model'], bp.P['gapfill'], bp.P['broken'])
bp.PARAM['light_uptake'] = 100.0; net = bp.autotrophic_bounds(m)
rows = list(csv.reader(open(os.path.join(bp.P['out'], f'biomreaction_{cond}.csv')))); hdr = rows[0]; lines = hdr[3:]
tab3 = {r['model_metabolite']: r for r in bp.read_csv(bp.P['bof'])}
# BOF vectors over model species (+ ATP maintenance), restricted to balanced rows
nS = m['S'].shape[0]; B = np.zeros((nS, len(lines)))
for r in rows[1:]:
    comp = r[1]; i = bp.species_for(m, r[0], (comp, 'c', 'p', 'm', 'x', 'v'))
    B[i] -= np.array(r[3:], float)
for nm, s in (('ATP', -1), ('ADP', 1), ('Orthophosphate', 1)):
    B[bp.species_for(m, nm, ('c',))] += s * bp.PARAM['atp_maintenance']
B = B[net['keep']]
S = np.hstack([net['S'], B[:, [lines.index('B73')]]]); lb = np.append(net['lb'], 0); ub = np.append(net['ub'], 1000); jb = S.shape[1] - 1
Z = bp.load_zhang(bp.P['zhang'], bp.P['cache'])
dm = np.array([Z['s1'].get(l, {}).get('dry') or np.nan for l in lines], float)
keep = np.where(np.isfinite(dm))[0][:NMAX]
t0 = time.time()
ck = os.environ.get('NETGS_CHUNK')   # e.g. "0:60" -> QP only for these line positions, saved to a part file
sub = None
if ck:
    a_, b_ = map(int, ck.split(':')); sub = set(range(a_, min(b_, len(keep))))
r = genotype_fluxes(S, lb, ub, jb, B[:, keep], dm[keep], ref=list(keep).index(lines.index('B73')) if lines.index('B73') in keep else 0,
                    objective=OBJ, ref_bof=B[:, lines.index('B73')], smodel=os.environ.get('NETGS_SMODEL', 'mean'), qp_subset=sub,
                    solver=SOLVER, extra_rxns=[net['rxns'].index(x) for x in EXTRA])
print(f'{len(keep)} lines, objective={OBJ}, solver={SOLVER}, {time.time()-t0:.0f}s; active reactions in reference: {len(r["nonzero"])}; '
      f's_model={r["s_model"]:.4g} (min {r["zmax"].min():.4g}); status: { {s: r["status"].count(s) for s in set(r["status"])} }')
rx = net['rxns'] + ['BIOMASS']
if ck:
    np.save(os.path.join(bp.P['out'], f'netgs_part_{cond}_{TAG}_{ck.replace(":","-")}.npy'), dict(fluxes=r['fluxes'], status=r['status'], nonzero=r['nonzero'], vref=r['vref'], target=r['target'], half=r['halfwidth'], s_model=r['s_model'], zmax=r['zmax'], keep=keep, lines=lines, rx=rx), allow_pickle=True)
    print('chunk', ck, 'done', time.time()-t0); sys.exit(0)
out = os.path.join(bp.P['out'], f'netgs_fluxes_{cond}_{TAG}.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['reaction', 'v_ref_B73'] + [lines[k] for k in keep])
    for a, j in enumerate(r['nonzero']):
        w.writerow([rx[j], f'{r["vref"][j]:.6g}'] + [f'{x:.6g}' for x in r['fluxes'][a]])
kb = int(np.where(r['nonzero'] == jb)[0][0])
print('biomass flux vs target: corr %.4f, max abs dev %.2e (window +/- %.2e)' % (np.corrcoef(r['fluxes'][kb], r['target'])[0, 1],
      np.nanmax(np.abs(r['fluxes'][kb] - r['target'])), r['halfwidth']))
json.dump(dict(objective=OBJ, n_lines=int(len(keep)), s_model=float(r['s_model']), zmax_min=float(r['zmax'].min()),
               halfwidth=float(r['halfwidth']), n_active=int(len(r['nonzero'])), status={s: r['status'].count(s) for s in set(r['status'])}),
          open(os.path.join(bp.P['out'], f'netgs_summary_{cond}_{TAG}.json'), 'w'), indent=1)
print('written', out)
