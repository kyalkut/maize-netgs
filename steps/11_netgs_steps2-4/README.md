# 11_netgs_steps2-4

Genotypes (ZEAMAP AMP SNPs, B73 v4) -> kernel K = ZZ'; netGS step 2 (rrBLUP per flux, 50 x 3-fold CV), step 3 (QP towards predicted fluxes with the test line's BOF), step 4 (r with DS dry mass) vs classical GS on dry mass; both BOF variants.

**Scripts:** snp_prepare.py; netgs_steps234.py (NETGS_TAG = step-1 run, NETGS_OBJ3 = code/paper, BOF_OUT = variant folder)

**Inputs:**
- `inputs/snp_line_match.csv` (from `output/snp_line_match.csv`)
- `inputs/netgs_foldid_DS.csv` (from `output/netgs_foldid_DS.csv`)

**Outputs:**
- `outputs/snp_summary.json` (from `output/snp_summary.json`)
- `outputs/snp_kernel_DS.npz` (from `output/snp_kernel_DS.npz`)
- `outputs/netgs_step2_accuracy_DS_paper_cvxpy_ser.csv` (from `output/netgs_step2_accuracy_DS_paper_cvxpy_ser.csv`)
- `outputs/netgs_step2_accuracy_DS_code_cvxpy_ser.csv` (from `output/netgs_step2_accuracy_DS_code_cvxpy_ser.csv`)
- `outputs/netgs_step2_accuracy_DS_paper_cvxpy_co.csv` (from `output/netgs_step2_accuracy_DS_paper_cvxpy_co.csv`)
- `outputs/netgs_step2_accuracy_DS_paper_cvxpy_ser.csv` (from `output/no24/netgs_step2_accuracy_DS_paper_cvxpy_ser.csv`)
- `outputs/netgs_step2_accuracy_DS_code_cvxpy_ser.csv` (from `output/no24/netgs_step2_accuracy_DS_code_cvxpy_ser.csv`)
- `outputs/netgs_step2_accuracy_DS_paper_cvxpy_co.csv` (from `output/no24/netgs_step2_accuracy_DS_paper_cvxpy_co.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_ser_paper.json` (from `output/netgs_step34_DS_paper_cvxpy_ser_paper.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_ser_paper.csv` (from `output/netgs_step34_folds_DS_paper_cvxpy_ser_paper.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_ser_code.json` (from `output/netgs_step34_DS_paper_cvxpy_ser_code.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_ser_code.csv` (from `output/netgs_step34_folds_DS_paper_cvxpy_ser_code.csv`)
- `outputs/netgs_step34_DS_code_cvxpy_ser_paper.json` (from `output/netgs_step34_DS_code_cvxpy_ser_paper.json`)
- `outputs/netgs_step34_folds_DS_code_cvxpy_ser_paper.csv` (from `output/netgs_step34_folds_DS_code_cvxpy_ser_paper.csv`)
- `outputs/netgs_step34_DS_code_cvxpy_ser_code.json` (from `output/netgs_step34_DS_code_cvxpy_ser_code.json`)
- `outputs/netgs_step34_folds_DS_code_cvxpy_ser_code.csv` (from `output/netgs_step34_folds_DS_code_cvxpy_ser_code.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_co_paper.json` (from `output/netgs_step34_DS_paper_cvxpy_co_paper.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_co_paper.csv` (from `output/netgs_step34_folds_DS_paper_cvxpy_co_paper.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_co_code.json` (from `output/netgs_step34_DS_paper_cvxpy_co_code.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_co_code.csv` (from `output/netgs_step34_folds_DS_paper_cvxpy_co_code.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_ser_paper.json` (from `output/no24/netgs_step34_DS_paper_cvxpy_ser_paper.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_ser_paper.csv` (from `output/no24/netgs_step34_folds_DS_paper_cvxpy_ser_paper.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_ser_code.json` (from `output/no24/netgs_step34_DS_paper_cvxpy_ser_code.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_ser_code.csv` (from `output/no24/netgs_step34_folds_DS_paper_cvxpy_ser_code.csv`)
- `outputs/netgs_step34_DS_code_cvxpy_ser_paper.json` (from `output/no24/netgs_step34_DS_code_cvxpy_ser_paper.json`)
- `outputs/netgs_step34_folds_DS_code_cvxpy_ser_paper.csv` (from `output/no24/netgs_step34_folds_DS_code_cvxpy_ser_paper.csv`)
- `outputs/netgs_step34_DS_code_cvxpy_ser_code.json` (from `output/no24/netgs_step34_DS_code_cvxpy_ser_code.json`)
- `outputs/netgs_step34_folds_DS_code_cvxpy_ser_code.csv` (from `output/no24/netgs_step34_folds_DS_code_cvxpy_ser_code.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_co_paper.json` (from `output/no24/netgs_step34_DS_paper_cvxpy_co_paper.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_co_paper.csv` (from `output/no24/netgs_step34_folds_DS_paper_cvxpy_co_paper.csv`)
- `outputs/netgs_step34_DS_paper_cvxpy_co_code.json` (from `output/no24/netgs_step34_DS_paper_cvxpy_co_code.json`)
- `outputs/netgs_step34_folds_DS_paper_cvxpy_co_code.csv` (from `output/no24/netgs_step34_folds_DS_paper_cvxpy_co_code.csv`)

The VCF (208 MB, input/snp/) is not copied: ftp.cngb.org/pub/CNSA/data3/CNP0001565/zeamap/02_Variants/SNPs/AMP_SNP_anno.vcf.gz.
