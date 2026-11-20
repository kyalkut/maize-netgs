# 10_variant_without_24

Variant: the 24 metabolites without a literature anchor left out (BOF = Table III + 7 amino acids, 37 components, 1 g DW); BOFs, netGS step 1 (paper/code x serine pathway/RuBisCO bounds), FluxSum, data analyses; compared with the main run (61 components).

**Scripts:** bof_pipeline.py with BOF_NO_PH=1 BOF_OUT=output/no24; run_netgs_maize.py + merge_netgs.py + fluxsum.py + analysis_data.py with BOF_OUT=output/no24; compare_variants.py

**Inputs:**
- `inputs/bof_components.csv` (from `output/no24/bof_components.csv`)
- `inputs/biomreaction_DS.csv` (from `output/no24/biomreaction_DS.csv`)

**Outputs:**
- `outputs/no24/` (from `output/no24`)
- `outputs/variant_comparison.json` (from `output/variant_comparison.json`)
- `outputs/variant_comparison_bof_DS.csv` (from `output/variant_comparison_bof_DS.csv`)
- `outputs/variant_comparison_reactions_paper_cvxpy_ser.csv` (from `output/variant_comparison_reactions_paper_cvxpy_ser.csv`)
- `outputs/variant_comparison_reactions_code_cvxpy_ser.csv` (from `output/variant_comparison_reactions_code_cvxpy_ser.csv`)
- `outputs/variant_comparison_reactions_paper_cvxpy_co.csv` (from `output/variant_comparison_reactions_paper_cvxpy_co.csv`)
- `outputs/variant_comparison_reactions_code_cvxpy_co.csv` (from `output/variant_comparison_reactions_code_cvxpy_co.csv`)

output/no24 holds the full set of files of the variant (same names as in output/).
