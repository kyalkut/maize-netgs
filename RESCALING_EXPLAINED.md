# Selection and rescaling: which metabolites enter the BOF, and how much

Two steps of one job: **A** chooses which measured metabolites get a data-driven (non-zero) coefficient; **B** computes how much of each, per line and condition. Running example: B73, tryptophan, drought (DS).

## Files in this folder, by step

| Step | File | What it holds |
|---|---|---|
| A1 | `metabolite_mapping.csv` | Zhang feature → C4GEM metabolite (52 rows); `include` = yes/no switches a match off |
| A2–A4 | `run_log.txt` (STEP 3 and 4) | Which metabolites were kept or dropped, and why |
| A result | `bof_components.csv` | All 61 recipe components; `varies = yes` marks the 31 data-driven ones; bound part, anchor, source, molecular weight |
| B1 | `relative_factors_WW.csv`, `relative_factors_DS.csv` | r per metabolite (rows) × line (columns) |
| B2 | `anchors.csv` | Absolute free-pool anchor per metabolite and its source (empty = 1 µmol/gDW placeholder) |
| B2 | `free_pool_mmol_WW.csv`, `free_pool_mmol_DS.csv` | Rescaled free amount = anchor × r, metabolite × line |
| B3–B4 | `biomreaction_WW.csv`, `biomreaction_DS.csv` | Final coefficients, 61 components × 385 lines (what netGS reads) |
| all of B | `rescaling_trace.csv` | Every step for every condition × line × metabolite (columns 1–9, below) |
| check | `qc_bof_coefficients.png` | Spread of each varying coefficient across lines; WW blue, DS orange |

## A. Which metabolites enter the recipe (351 → 31)

1. **In the model?** Name matching of Zhang's 351 annotated features to C4GEM (normalised names + ~20 manual synonyms, look-alikes rejected): **52**.
2. **Producible?** The model must make each component from CO₂, light, N, S, P: **35** (no route for raffinose, riboflavin, choline, pantothenate, NAD by-products, …).
3. **Clean data?** Measured in ≥ 90% of lines in both conditions; SD of log2 ≤ 2.5; one feature per metabolite (most complete, then most intense): **31** (Arg absent under DS; GMP/UMP incomplete; AMP noisy; duplicate 4-coumarate).
4. **Cleaning:** clip to the 1st–99th percentile; fill gaps with the condition mean.

Result: 7 Table III amino acids (Glu, Leu, Met, Phe, Thr, Trp, Tyr) + 24 new (indole, citrate, isocitrate, nucleosides, ferulate, …) vary by line; 30 unmeasured Table III components stay fixed. 61 components in total.

## B. How much (rescaling)

| | Formula | B73, Trp, DS |
|---|---|---|
| B1 ratio | r = 2^(log2 value − WW mean over 385 lines) | 2^(17.56 − 14.05) = **11.4** |
| B2 free part | anchor × r | 0.000435 × 11.4 = **0.00495 mmol** |
| B3 two pools | bound (Table III, fixed; 0 if not in Table III) + free | 0.040 + 0.00495 = 0.045 |
| B4 equal mass | fixed part × s, s = (536.3 − free mass) / fixed mass | 0.040 × 0.993 + 0.00495 = **0.0447** |

- **WW mean for both conditions** so the drought shift is kept (Trp median r: 1.0 WW, 4.5 DS).
- **Anchors:** 7 amino acids from Shu et al. 2022 (maize leaf, V6, LC-MS; µmol/L × 0.05 = µmol/gDW, assuming 20 mg dry powder in 1 mL); 24 placeholders at 1 µmol/gDW.
- **Two pools** because Table III amino acids are mostly protein-bound while Zhang measures free pools; scaling the whole coefficient would give 11× Trp in protein.
- **Equal mass:** target = average recipe, 536.3 mg (Table III alone = 529 mg, about half of dry mass); s between 0.98 and 1.00; free parts untouched.

## `rescaling_trace.csv` columns

| # | Column | Meaning |
|---|---|---|
| 1 | `1_log2_intensity` | Zhang Table S3 value (after cleaning) |
| 2 | `2_mean_log2_WW` | WW mean over lines (reference) |
| 3 | `3_relative_factor` | r = 2^(1 − 2) |
| 4 | `4_free_anchor_mmol` | anchor (from `anchors.csv` or placeholder) |
| 5 | `5_free_part_mmol` | 4 × 3: the rescaled value |
| 6 | `6_bound_part_tableIII` | fixed protein-bound part |
| 7 | `7_fixed_part_scale` | s, mass equalisation |
| 8 | `8_bound_part_scaled` | 6 × 7 |
| 9 | `9_final_coefficient` | 5 + 8, identical to `biomreaction_*.csv` |

## Checks

- All 385 lines' recipes are feasible in both conditions.
- Placeholder anchors ×0.1 / ×10: spread between lines 0.2% → 7%, line ranking r = 0.74–0.98. Anchors set the size of differences, not their direction.

## Choices open for discussion

Name-based matching; 90% coverage and SD ≤ 2.5 thresholds; WW mean as reference; unit conversion of Shu et al.; mass equalised to the average recipe rather than 1 g; 24 placeholder anchors. All are set in `curation/` or in `PARAM` in `bof_pipeline.py`.
