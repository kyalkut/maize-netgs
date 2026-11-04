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
# RuBisCO carboxylation/oxygenation bounds lo <= vc/vo <= hi on the plastid pair, e.g. "22.2:88.0" (maize, C4GEM;
# Wang et al. 2012 BMC Syst Biol 6 S2:S9, 100-1000 ubar CO2); netGS Arabidopsis used 0.938:3.814 ('none' = off)
CO = os.environ.get('NETGS_CO', 'none')
if CO.lower() != 'none': TAG += '_co'
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
extra_ub = None
if CO.lower() != 'none':
    lo_r, hi_r = map(float, CO.split(':')); jc, jo = net['rxns'].index('R_R00024_p'), net['rxns'].index('R_R03140_p')
    A = np.zeros((2, S.shape[1])); A[0, jc], A[0, jo] = 1, -hi_r; A[1, jc], A[1, jo] = -1, lo_r   # vc - hi*vo <= 0, lo*vo - vc <= 0
    extra_ub = (A, np.zeros(2))
t0 = time.time()
ck =os.environ.get('NETGS_CHUNK')   # e.g. "0:60" -> QP only for these line positions, saved to a part file
sub = None
if ck:
    a_, b_ = map(int, ck.split(':')); sub = set(range(a_, min(b_, len(keep))))
r = genotype_fluxes(S, lb, ub, jb, B[:, keep], dm[keep], ref=list(keep).index(lines.index('B73')) if lines.index('B73') in keep else 0,
                    objective=OBJ, ref_bof=B[:, lines.index('B73')], smodel=os.environ.get('NETGS_SMODEL', 'mean'), qp_subset=sub,
                    solver=SOLVER, extra_rxns=[net['rxns'].index(x) for x in EXTRA], extra_ub=extra_ub)
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
