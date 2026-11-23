"""Genotypes for netGS steps 2-4: SNPs of the maize association panel (AMP) -> genomic relationship kernel of the
Zhang et al. 2021 lines with DS dry mass.

Source: ZEAMAP (Gui et al. 2020, iScience 23:101241), AMP SNPs from ~20x resequencing of 507 lines, B73 AGPv4,
  ftp.cngb.org/pub/CNSA/data3/CNP0001565/zeamap/02_Variants/SNPs/AMP_SNP_anno.vcf.gz  (md5 b95c23d3f2f788e7e1550743f59a7a28)
  -> input/snp/ (not in git). maizego.org lists the 1.25 M / 540-line set only on Baidu Pan (needs an account).

Line names: Zhang's names upper-cased, non-alphanumerics removed; if absent, the same name with suffix 'X' (the AMP's
  re-sampled entries, e.g. CIMBL45X; listed in snp_line_match.csv as 'x_suffix').
SNPs: biallelic, call rate >= SNP_CALL (0.8) and MAF >= SNP_MAF (0.05) within the matched lines.
Coding as rrBLUP: 0/0 -> -1, 0/1 -> 0, 1/1 -> 1; missing -> SNP mean.
rrBLUP's mixed.solve(y, Z = SNPs, K = NULL) assumes u ~ N(0, I s2u) per marker, i.e. genetic values g = Z u ~ N(0, s2u ZZ').
  So only K = ZZ' is needed (identical predictions, GBLUP form); it is accumulated chunk by chunk and saved:
  output/snp_kernel_DS.npz (K, lines, n_snps), output/snp_line_match.csv, output/snp_summary.json."""
import os, sys, gzip, re, csv, json, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
VCF = os.environ.get('SNP_VCF', os.path.join(HERE, 'input', 'snp', 'AMP_SNP_anno.vcf.gz'))
OUT = os.path.join(HERE, 'output')
MAF, CALL = float(os.environ.get('SNP_MAF', '0.05')), float(os.environ.get('SNP_CALL', '0.8'))
norm = lambda s: re.sub(r'[^A-Z0-9]', '', s.upper())

feas = list(csv.DictReader(open(os.path.join(OUT, 'feasibility_per_line.csv'))))
zl = [r['line'] for r in feas if r['dry_mass_DS_g']]          # lines with DS dry mass (netGS step 1 lines)
t0 = time.time()
with gzip.open(VCF, 'rt') as f:
    for line in f:
        if line.startswith('#CHROM'): samples = line.rstrip('\n').split('\t')[9:]; break
vn = {norm(s): i for i, s in enumerate(samples)}
match, how = {}, {}
for l in zl:
    k = norm(l)
    if k in vn: match[l] = vn[k]; how[l] = 'exact'
    elif k + 'X' in vn: match[l] = vn[k + 'X']; how[l] = 'x_suffix'
lines = [l for l in zl if l in match]; col = np.array([match[l] for l in lines]) + 9
print(f'{len(samples)} VCF samples; {len(zl)} Zhang lines with DS dry mass; matched {len(lines)} '
      f'({sum(how[l] == "x_suffix" for l in lines)} via X suffix)', flush=True)
with open(os.path.join(OUT, 'snp_line_match.csv'), 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['line', 'vcf_sample', 'match'])
    for l in zl: w.writerow([l, samples[match[l]] if l in match else '', how.get(l, 'not in AMP VCF')])

code = {'0/0': -1.0, '0|0': -1.0, '0/1': 0.0, '0|1': 0.0, '1|0': 0.0, '1/1': 1.0, '1|1': 1.0}
n = len(lines); K = np.zeros((n, n)); buf = []; n_tot = n_kept = 0; per_chr = {}
def flush():
    global K, buf
    if buf: Z = np.array(buf, dtype=np.float32).T; K += (Z.astype(np.float64) @ Z.T.astype(np.float64)); buf = []
with gzip.open(VCF, 'rt') as f:
    for line in f:
        if line[0] == '#': continue
        n_tot += 1
        p = line.split('\t', 9 + 1)                            # cheap check of REF/ALT first
        if len(p[3]) != 1 or len(p[4]) != 1: continue
        p = line.rstrip('\n').split('\t')
        g = np.array([code.get(p[c].split(':', 1)[0], np.nan) for c in col])
        ok = np.isfinite(g); cr = ok.mean()
        if cr < CALL: continue
        af = (g[ok].mean() + 1) / 2; maf = min(af, 1 - af)
        if maf < MAF: continue
        g[~ok] = g[ok].mean(); buf.append(g); n_kept += 1; per_chr[p[0]] = per_chr.get(p[0], 0) + 1
        if len(buf) >= 20000: flush()
        if n_tot % 1000000 == 0: print(f'  {n_tot} variants read, {n_kept} kept, {time.time()-t0:.0f}s', flush=True)
flush()
np.savez_compressed(os.path.join(OUT, 'snp_kernel_DS.npz'), K=K, lines=np.array(lines), n_snps=n_kept)
S = dict(vcf=os.path.basename(VCF), source='ZEAMAP AMP SNPs (Gui et al. 2020), B73 AGPv4', n_vcf_samples=len(samples),
         n_zhang_lines_with_drymass=len(zl), n_matched=len(lines), n_matched_x_suffix=sum(how[l] == 'x_suffix' for l in lines),
         not_matched=[l for l in zl if l not in match], n_variants_read=n_tot, n_snps_kept=n_kept, maf_min=MAF, call_rate_min=CALL,
         snps_per_chromosome=per_chr, coding='-1/0/1 (rrBLUP), missing = SNP mean', kernel='K = Z Z^T (equivalent to rrBLUP mixed.solve with Z, K = NULL)')
json.dump(S, open(os.path.join(OUT, 'snp_summary.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in S.items() if k != 'not_matched'}, indent=1), f'{time.time()-t0:.0f}s')
