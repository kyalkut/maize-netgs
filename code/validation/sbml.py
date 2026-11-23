"""Minimal SBML (COBRA level-2 style) reader -> S, lb, ub, c, ids. No external deps."""
import re, html, numpy as np
def read_model(path):
    t = open(path, encoding='utf-8').read()
    sp = re.findall(r'<species ([^>]*)>', t)
    mets = [re.search(r'id="([^"]+)"', s).group(1) for s in sp]
    names = [html.unescape((re.search(r'name="([^"]*)"', s) or re.search(r'id="([^"]+)"', s)).group(1)) for s in sp]
    midx = {m: i for i, m in enumerate(mets)}
    rxs = re.findall(r'<reaction ([^>]*)>(.*?)</reaction>', t, re.S)
    n = len(rxs); S = np.zeros((len(mets), n)); lb = np.zeros(n); ub = np.zeros(n); c = np.zeros(n); rids = []; rnames = []
    for j, (attr, body) in enumerate(rxs):
        rids.append(re.search(r'id="([^"]+)"', attr).group(1))
        nm = re.search(r'name="([^"]*)"', attr); rnames.append(html.unescape(nm.group(1)) if nm else rids[-1])
        rev = 'reversible="false"' not in attr
        for tag, sign in (('listOfReactants', -1), ('listOfProducts', 1)):
            m = re.search(rf'<{tag}>(.*?)</{tag}>', body, re.S)
            if not m: continue
            for ref in re.findall(r'<speciesReference ([^>]*?)/?>', m.group(1)):
                s = re.search(r'species="([^"]+)"', ref).group(1)
                st = re.search(r'stoichiometry="([^"]+)"', ref)
                S[midx[s], j] += sign * (float(st.group(1)) if st else 1.0)
        par = dict(re.findall(r'<parameter id="([A-Z_]+)" value="([^"]+)"', body))
        lb[j] = float(par.get('LOWER_BOUND', -1000 if rev else 0)); ub[j] = float(par.get('UPPER_BOUND', 1000))
        c[j] = float(par.get('OBJECTIVE_COEFFICIENT', 0))
    return dict(S=S, lb=lb, ub=ub, c=c, mets=mets, metNames=names, rxns=rids, rxnNames=rnames)
