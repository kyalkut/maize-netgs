"""Cross-check of maize models: C4GEM (Dal'Molin 2010, used here) vs Simons et al. 2014 (Plant Physiol 166:1659, two-cell
leaf model) vs Bogart & Myers 2016 (PLoS One e0151722; iEB5204 single cell, iEB2140x2 two-cell, from CornCyc).
Size, compartments, genes (GPR coverage, ID type), biomass reaction, overlap with Zhang's 351 annotated metabolites
(the same automatic name matching for every model). Model files: input/other_models/ (see README there)."""
import os, re, sys, csv, json, html, collections
import xml.etree.ElementTree as ET
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bof_pipeline as bp
OM = os.path.join(HERE, 'input', 'other_models'); OUT = bp.P['out']
MODELS = {'C4GEM': bp.P['model'],
          'Simons2014': os.path.join(OM, 'simons2014', 'sbml', 'TS4_Model.xml'),
          'iEB5204': os.path.join(OM, 'pone.0151722.s019.xml'),
          'iEB2140x2': os.path.join(OM, 'pone.0151722.s021.xml')}

def tag(e): return e.tag.rsplit('}', 1)[-1]
def read(path):
    root = ET.parse(path).getroot(); model = next(e for e in root if tag(e) == 'model')
    sp, rx = {}, []
    for e in model.iter():
        t = tag(e)
        if t == 'species':
            sp[e.get('id')] = dict(name=html.unescape(e.get('name') or e.get('id')), comp=e.get('compartment'),
                                   boundary=e.get('boundaryCondition') == 'true')
        elif t == 'reaction':
            st = collections.defaultdict(float); gpr = ''; lb = ub = None; obj = 0.0
            for ch in e.iter():
                c = tag(ch)
                if c in ('listOfReactants', 'listOfProducts'):
                    sg = -1 if c == 'listOfReactants' else 1
                    for s in ch:
                        st[s.get('species')] += sg * float(s.get('stoichiometry') or 1)
                elif c == 'p' and ch.text and ch.text.strip().startswith('GENE_ASSOCIATION:'):
                    gpr = ch.text.split(':', 1)[1].strip()
                elif c == 'p' and ch.text and ch.text.strip().startswith('Gene Id:'):            # C4GEM: one line per gene
                    g = ch.text.split(':', 1)[1].strip()
                    if g.isdigit(): gpr = (gpr + ' or ' if gpr else '') + g
                elif c == 'parameter':
                    pid, v = ch.get('id'), float(ch.get('value'))
                    if pid == 'LOWER_BOUND': lb = v
                    elif pid == 'UPPER_BOUND': ub = v
                    elif pid == 'OBJECTIVE_COEFFICIENT': obj = v
            rev = e.get('reversible') != 'false'
            rx.append(dict(id=e.get('id'), name=html.unescape(e.get('name') or ''), st=dict(st), gpr=gpr,
                           lb=lb if lb is not None else (-1000.0 if rev else 0.0), ub=ub if ub is not None else 1000.0,
                           obj=obj, bounds_in_file=lb is not None))
    return sp, rx

def norm(s):
    s = html.unescape(s or '').lower()
    s = re.sub(r'^(ms|bs)_', '', s)
    s = re.sub(r'\s*\[[^\]]*\]$', '', s)                       # "[biomass]" suffixes
    s = re.sub(r'\s*\([a-z0-9_+\-]+\)$', '', s) if '(' in s and not s.startswith('(') else s   # CornCyc "(CPD-ID)" suffix
    s = re.sub(r'_(c|m|p|v|x|e|ext|acc|pmd|total|cytosol|chloroplast|mitochondrion|peroxisome|vacuole)$', '', s)
    s = re.sub(r'^(l|d|s|r)-(?=[a-z])', '', s)                 # stereo prefixes
    s = re.sub(r"[\s\-_,'()]+", '', s)
    return s.replace('acid', 'ate')

def genes(gpr): return set(g for g in re.findall(r'[A-Za-z0-9_.\-]+', gpr) if g.lower() not in ('and', 'or'))

Z = bp.load_zhang(bp.P['zhang'], bp.P['cache'])
ann = {fid: d['name'] for fid, d in Z['s2'].items() if d.get('name') and not d['name'].startswith(('Unknown', 'unknown'))}
annorm = collections.defaultdict(set)
for fid, n in ann.items(): annorm[norm(n)].add(fid)
bof_ids = {c['zhang_feature'] for c in csv.DictReader(open(os.path.join(OUT, 'bof_components.csv'))) if c['varies'] == 'yes'}

res = {}
for name, path in MODELS.items():
    sp, rx = read(path)
    names = {norm(s['name']) for s in sp.values()}
    comps = collections.Counter(s['comp'] for s in sp.values())
    cnames = ' '.join(e.get('name') or '' for e in ET.parse(path).getroot().iter() if tag(e) == 'compartment')
    cells = sorted({p for p in ('ms_', 'bs_') if any(k.startswith(p) for k in sp)}) or             (['bundle sheath', 'mesophyll'] if re.search('bundle', cnames, re.I) and re.search('mesophyll', cnames, re.I) else ['single'])
    gs = set().union(*[genes(r['gpr']) for r in rx]) if rx else set()
    gs_base = {re.sub(r'^(ms|bs)_', '', g) for g in gs}
    idtype = collections.Counter('GRMZM/AC (B73 v2-v3)' if re.match(r'(GRMZM|AC\d)', g) else 'Zm00001 (v4+)' if g.startswith('Zm0000')
                                 else 'numeric (NCBI Entrez)' if g.isdigit() else 'other' for g in gs_base)
    bio = [r for r in rx if re.search(r'biomass', r['name'] + r['id'], re.I) and len(r['st']) > 5]
    biomain = max(bio, key=lambda r: len(r['st'])) if bio else None
    if name == 'C4GEM':      # BOF only in Dal'Molin Table III, not in the SBML
        bio_comps = sorted({norm(r['model_metabolite']) for r in bp.read_csv(bp.P['bof'])})
    elif any('biomass' in (s['comp'] or '').lower() for s in sp.values()):   # Bogart: components live in a biomass compartment
        bio_comps = sorted({norm(s['name']) for s in sp.values() if 'biomass' in (s['comp'] or '').lower()})
    else:
        bio_comps = sorted({norm(sp[s]['name']) for s, v in biomain['st'].items() if v < 0}) if biomain else []
    hits = {k for k in annorm if k in names}; feats = set().union(*[annorm[k] for k in hits]) if hits else set()
    res[name] = dict(file=os.path.relpath(path, HERE), reactions=len(rx), species=len(sp), unique_metabolite_names=len(names),
                     compartments=len(comps), cell_types=cells, reactions_with_gpr=sum(1 for r in rx if r['gpr']),
                     genes=len(gs_base), gene_id_types=dict(idtype), bounds_in_file=any(r['bounds_in_file'] for r in rx),
                     biomass_reactions=[(r['id'], len(r['st'])) for r in bio][:6],
                     biomass_main='Table III (paper)' if name == 'C4GEM' else (None if biomain is None else biomain['id']), biomass_components=len(bio_comps),
                     zhang_annotated_matched=len(feats), zhang_names_matched=len(hits),
                     zhang_bof31_matched=len(feats & bof_ids),
                     zhang_in_biomass=sorted(k for k in hits if k in set(bio_comps)),
                     _names=names, _hits=hits)
    print(f"{name:10} rxns {len(rx):5}  species {len(sp):5}  comps {len(comps):2}  cells {cells}  GPR rxns {res[name]['reactions_with_gpr']:5}  "
          f"genes {len(gs_base):5} {dict(idtype)}  biomass {res[name]['biomass_main']} ({len(bio_comps)} comps)  "
          f"Zhang annotated matched {len(feats)}  of BOF-31 {len(feats & bof_ids)}")
# which of our 31 BOF metabolites are in each model; union/intersection of Zhang matches
allhits = {k: v['_hits'] for k, v in res.items()}
rows = []
for key in sorted(set().union(*allhits.values())):
    fids = sorted(annorm[key]); rows.append(dict(zhang_name=ann[fids[0]], zhang_features=';'.join(fids), in_our_bof='yes' if set(fids) & bof_ids else 'no',
                                                 **{m: 'yes' if key in allhits[m] else '' for m in MODELS}))
with open(os.path.join(OUT, 'model_comparison_zhang_overlap.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
json.dump({k: {a: b for a, b in v.items() if not a.startswith('_')} for k, v in res.items()},
          open(os.path.join(OUT, 'model_comparison_summary.json'), 'w'), indent=1)
print('Zhang names matched by any model:', len(rows), '| only Simons/Bogart, not C4GEM:',
      sum(1 for r in rows if not r['C4GEM'] and (r['Simons2014'] or r['iEB5204'] or r['iEB2140x2'])))
# FBA check of the Simons SBML as published (file bounds, light <= 100 per cell): does growth need CO2?
import numpy as np, scipy.sparse as sps
from scipy.optimize import linprog
spc, rx = read(MODELS['Simons2014']); ids = [k for k, v in spc.items() if not v['boundary']]; ix = {k: i for i, k in enumerate(ids)}
tr = [(ix[s], j, v) for j, r in enumerate(rx) for s, v in r['st'].items() if s in ix]
S = sps.csr_matrix(([t[2] for t in tr], ([t[0] for t in tr], [t[1] for t in tr])), shape=(len(ids), len(rx)))
rid = [r['id'] for r in rx]; lb = np.array([r['lb'] for r in rx]); ub = np.array([r['ub'] for r in rx])
for k in ('R_ExBe15', 'R_ExMe15'): lb[rid.index(k)], ub[rid.index(k)] = -100, 100
c = np.zeros(len(rx)); c[rid.index('R_Bio_Nplus')] = -1
sol = linprog(c, A_eq=S, b_eq=np.zeros(S.shape[0]), bounds=list(zip(lb, ub)), method='highs')
fx = {k: float(sol.x[rid.index(k)]) for k in ('R_Bio_Nplus', 'R_ExMe4', 'R_ExMe17', 'R_ExBe8', 'R_R00024_p', 'R_R03140_p', 'R_R00345_d', 'R_R00206_q')}
print('Simons FBA (file bounds, light 100/cell):', {k: round(v, 3) for k, v in fx.items()})
s = json.load(open(os.path.join(OUT, 'model_comparison_summary.json')))
s['Simons2014']['fba_file_bounds_light100'] = dict(fluxes=fx, note='growth without CO2 uptake (R_ExMe4/R_ExMe17 = 0) while maltose is exported (R_ExBe8): '
                                                   'the published SBML bounds are not the paper simulation setup; no fair FBA comparison from these files')
s['iEB5204']['fba_note'] = s['iEB2140x2']['fba_note'] = 'no flux bounds in the SBML; constraints are in S2 Protocol (fluxtools), and the leaf model is nonlinear'
json.dump(s, open(os.path.join(OUT, 'model_comparison_summary.json'), 'w'), indent=1)
