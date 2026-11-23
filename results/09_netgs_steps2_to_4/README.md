# 9. netGS steps 2 to 4

Scripts (in `code/`): `snp_prepare.py`, `netgs_steps234.py`, `analysis_data.py`

- `snp_summary.json`, `snp_line_match.csv`: markers kept (796,846) and how the 353 lines were matched (`match`).
  The marker file itself (208 MB) is not copied: ZEAMAP, ftp.cngb.org/pub/CNSA/data3/CNP0001565/zeamap/02_Variants/SNPs/AMP_SNP_anno.vcf.gz.
- `netgs_foldid_DS.csv`: the 50 x 3-fold split, the same for every run.
- `netgs_step2_accuracy_DS_paper_cvxpy_ser.csv`: per flux `accuracy_mean` (step 2), `r_with_biomass_flux`.
- `netgs_step34_DS_<step 1 run>_<step 3 objective>.json`: `netgs_accuracy_mean`, `classical_gs_drymass_mean` (0.325),
  `share_folds_netgs_better`. Suffix `_sel-<mode>`: which reactions enter step 3 (biomass, growth, all, inner, nongrowth).
  `no24/`: without the 24 placeholders. `sensitivity_x*/`: placeholders 10 times smaller or larger.
- `netgs_step34_folds_*.csv`: the same per test set (`r_netgs`, `r_classical_gs_drymass`).
- `analysis_summary.json`, `prediction`: metabolites predicting dry mass (all 0.54, 31 random 0.29, our 31 0.19).
