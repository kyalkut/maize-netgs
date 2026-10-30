"""
Curation of the published C4GEM maize SBML (Dal'Molin et al. 2010, Supplemental File S1).

Problems found in the published file and how they are fixed:
  A. Gap-filling reactions listed in the paper (Table S1) are missing from the SBML (57 of 100).
     -> 43 with parsable equations are added (+3 manual links: H2O2 c->x, UMP kinase, cytosolic fumarase),
        irreversible in the written direction unless marked reversible  (curation/gapfill_reactions.json)
  B. 12 reactions create or destroy carbon compounds from water because the generic polymer/peptide
     substrate was stripped when the model was exported (e.g. 'H2O -> D-Glucose', 'H2O -> L-Proline').
     -> blocked                                                         (curation/broken_reactions.csv)
  C. The SBML has no flux bounds; every reaction without reversible="false" is reversible, which creates
     energy-generating cycles (ATP made without light).
     -> standard directionality rules by EC class (kinases: ATP-consuming direction; phosphatases /
        anhydride hydrolases / amidases / glycosidases: hydrolysis direction; aldehyde dehydrogenases:
        oxidation direction), with exceptions for truly reversible enzymes (adenylate kinase, NDP kinase,
        phosphoglycerate kinase, GAPDH, ...), plus ATP hydrolysis R00086 irreversible.
After curation: no ATP or biomass without light, no biomass without CO2, all 37 Table III components producible.
"""
import csv, json, re
import numpy as np

REV_EXCEPT_EC = {'2.7.4.3', '2.7.4.6', '2.7.2.3', '2.7.4.8', '2.7.4.22', '2.7.4.14', '2.7.4.4', '2.7.4.9',
                 '1.2.1.12', '1.2.1.13', '1.2.1.59', '1.2.1.9', '1.2.1.11', '1.2.1.38', '1.2.1.41'}

def ec_numbers(sbml_path):
    t = open(sbml_path, encoding='utf-8').read(); d = {}
    for rid, body in re.findall(r'<reaction id="([^"]+)"[^>]*>(.*?)</reaction>', t, re.S):
        e = re.search(r'EC: *([0-9.\-; ,]+)', body)
        d[rid] = [x.strip() for x in re.split(r'[;, ]+', e.group(1)) if x.strip()] if e else []
    return d

def add_gapfill(model, path):
    gf = json.load(open(path)); idx = {s: i for i, s in enumerate(model['mets'])}; cols = []
    for r in gf:
        c = np.zeros(len(model['mets']))
        for co, s in r['L']: c[idx[s]] -= co
        for co, s in r['R']: c[idx[s]] += co
        cols.append(c)
    model['S'] = np.hstack([model['S'], np.array(cols).T])
    model['lb'] = np.append(model['lb'], [-1000.0 if r.get('reversible') else 0.0 for r in gf])
    model['ub'] = np.append(model['ub'], [1000.0] * len(gf))
    model['rxns'] = model['rxns'] + ['GF_' + r['id'] for r in gf]
    return len(gf)

def block_broken(model, path):
    rows = list(csv.DictReader(open(path, encoding='utf-8')))
    for r in rows:
        j = model['rxns'].index(r['reaction']); model['lb'][j] = 0; model['ub'][j] = 0
    return rows

def direction_rules(model, ecs):
    S = model['S']
    sets = {k: {i for i, b in enumerate(model['base']) if b in v} for k, v in
            dict(ATP={'ATP'}, H2O={'H2O'}, NADox={'NAD+', 'NADP+'}).items()}
    changed = []
    j = model['rxns'].index('R_R00086_c'); model['lb'][j] = 0; changed.append(('R_R00086_c', 'ATP hydrolysis: irreversible'))
    for j, r in enumerate(model['rxns']):
        if r not in ecs or not (model['lb'][j] < 0 < model['ub'][j]): continue
        E = ecs[r]
        if not E or any(e in REV_EXCEPT_EC for e in E): continue
        col = S[:, j]
        side = lambda s: 'L' if any(col[i] < 0 for i in s) else ('R' if any(col[i] > 0 for i in s) else None)
        rule = None
        if any(e.startswith(('2.7.1.', '2.7.2.', '2.7.4.', '2.7.7.')) for e in E) and side(sets['ATP']):
            rule = ('kinase: ATP-consuming direction', side(sets['ATP']))
        elif any(e.startswith(('3.1.3.', '3.6.1.', '3.1.1.', '3.5.1.', '3.2.1.')) for e in E) and side(sets['H2O']):
            rule = ('hydrolase: hydrolysis direction', side(sets['H2O']))
        elif any(e.startswith('1.2.1.') for e in E) and side(sets['NADox']):
            rule = ('aldehyde dehydrogenase: oxidation direction', side(sets['NADox']))
        if rule:
            if rule[1] == 'L': model['lb'][j] = 0
            else: model['ub'][j] = 0
            changed.append((r, rule[0]))
    return changed

def curate(model, sbml_path, gapfill_path, broken_path):
    n_gf = add_gapfill(model, gapfill_path)
    broken = block_broken(model, broken_path)
    rules = direction_rules(model, ec_numbers(sbml_path))
    return dict(gapfill=n_gf, broken=len(broken), direction_changes=rules)
