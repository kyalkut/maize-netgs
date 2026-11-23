# netGS in maize: progress since last meeting

**Your requests from the last meeting and where they stand**

| Your note | Status |
|---|---|
| Non-zero coefficients from the data | **Done.** The BOF is built from the data: 31 measured metabolites vary by line (previously 7 Table III matches) |
| Is the maize model specific to one line? | **Answered.** C4GEM is built from the B73 genome, so B73 is used as the reference line |
| C4GEM maize gene mapping: multiple / lacking (CHECK) | **Checked.** 648/1,243 reactions have genes, 595 have none; up to 62 genes per reaction; IDs are 2009 NCBI Entrez and must be converted to B73 v3 (`GRMZM…`) for SNPs |
| Biomass rescaling (!!!) / update netGS with the rescaling | **Done:** two-pool rescaling of relative data to BOF coefficients, plus growth rescaling (s_model, s_measurement) in a Python port of netGS |
| Update the pipeline | **Done:** runnable pipeline (data → BOFs → netGS step 1), validated on the original Arabidopsis example |
| One set of data | Not changed. Zhang only; pipeline produces WW and DS. *Please confirm what was meant* |

---

## Findings

### 1. Rescaling relative metabolomics into BOF coefficients (two-pool)
- Zhang's data are relative (log2 intensity) and measure *free* pools; BOF coefficients are absolute and mostly *protein-bound*.
- Scaling whole coefficients would give impossible protein compositions (Trp ×12, Met ×0.03 in some lines).
- **Two-pool BOF:** `b = b_bound (Table III, fixed) + b_free (literature anchor) × 2^(x_line − mean_WW)`.
  - Common WW anchor for both conditions, so the drought shift is kept.
  - Renormalised so every line has the same BOF mass.
- **Anchors:** 7 amino acids from maize leaf LC-MS data (Shu et al. 2022, inbred lines at V6, the same stage as Zhang). The other 24 are placeholders.
- **Sensitivity:** with placeholders ×0.1 / ×1 / ×10, every line stays feasible. The spread between lines scales with the anchors (0.2% / 0.8% / 7%), but **the ranking of lines is stable** (rank correlation 0.74–0.98).

### 2. The published C4GEM file is not a working model (new)
- **57 of the 100 gap-filling reactions** listed in the paper's own Table S1 are **missing from the SBML**. As published, the model cannot make any nucleotide (ATP, GTP, …). Added: 43 parsable ones + 3 missing links (H₂O₂ transport, UMP kinase, cytosolic fumarase).
- **12 reactions create or destroy compounds from water**: generic substrates were stripped on export, e.g. `H2O → D-Glucose`, `H2O → D-Xylose`, `H2O → L-Proline`. Blocked.
- **No flux bounds; most reactions reversible → energy loops**: the model makes ATP without light. Fixed with standard EC-based direction rules (kinases, hydrolases, aldehyde dehydrogenases; 28 reactions + ATP hydrolysis).
- **After curation:** no ATP or biomass without light, no biomass without CO₂, biomass scales with light, all 37 Table III components producible.
- **Table III covers only ~53% of dry mass** (529 mg per unit), so BOF mass is normalised to the reference BOF, not to 1 g.

### 3. Which measured metabolites can enter the BOF
- 52 matched features → 35 producible by the curated model → **31 kept** after data QC.
- Not producible (no synthesis route in C4GEM): raffinose, quercetin glucoside, riboflavin, pantothenate, choline, betaine, NAD-salvage intermediates, …
- Dropped for data reasons: Arg (not detected under DS), AMP/GMP/UMP (too noisy or incomplete), 2 duplicate 4-coumarate features.

### 4. Growth rescaling ("two values") is in place
- **s_model** = max biomass each line's BOF allows (curated model, light-limited); **s_measurement** = max DS dry mass. Target = (s_model − δ) × dry mass / max dry mass, window ±1.645 SEM.
- With the curated model, all 371 lines with growth data have reachable targets with either the mean (netGS default) or the minimum of the line maxima. I used the minimum to be safe.
- **Small signal:** the model's maximum biomass from the BOF alone correlates weakly with measured DS dry mass (r ≈ 0.10, p ≈ 0.05, same at all anchor scales).

### 5. netGS code vs. paper: the objective differs (please check)
- **Paper, Methods Eq. 2:** minimise Σᵢ [(v_i^Col0 − v_i^Z) / v_i^Col0]², the per-reaction relative squared deviation from the reference.
- **Published `1_FluxDist.m`:** `A = 1./fluxc.'` is a 1×n row vector, so `norm(A*x − b)` is a scalar, |Σᵢ v_i/v_i^Col0 − 1|.
- **Port validated on the authors' Arabidopsis example:** s_model = 0.0035036, identical to the value hard-coded in `1_FluxDist.m`.
- **The two objectives give different flux traits.** Per-reaction correlation across the 67 accessions, paper vs. code version: **median r = 0.25**. The code version is deterministic (two LP algorithms give the same answer), so it is not solver noise.
- Since netGS step 2 predicts exactly these per-reaction fluxes, this matters for the results. **Which version was used for the published accuracies?**
- (Minor: my FBA reference gives 319 active reactions vs. 336 in the code. Different solver, same optimum; FBA optima are not unique.)

### 6. G×E signal in the BOF metabolites (new)
- Line ranking WW vs DS: **BOF metabolites median r = 0.39 vs. 0.63 for all 3,116 metabolites.**
- The **free amino acids almost completely reorder** under drought (Glu 0.05, Met 0.06, Leu 0.07, Thr 0.10, Tyr 0.12). Strong G×E/plasticity candidates for objective 3. *Caveat: without replicates, G×E can't be separated from measurement noise.*

### 7. Growth information sits mostly outside the model (new, important for expectations)
Prediction of DS dry mass from metabolites (ridge, 5-fold CV, 371 lines):

| Metabolite set | r |
|---|---|
| 31 BOF metabolites | **0.19** |
| 31 random metabolites | 0.29 (95% range 0.19–0.36) |
| 825 drought-responsive | 0.52 |
| all 3,134 | 0.54 |

- The metabolites the model can use carry **less** growth information than random ones. The predictive signal is mostly in secondary metabolism, which C4GEM doesn't cover.
- Implication: netGS in maize may gain less from line-specific BOFs than in Arabidopsis. A fair comparison with classical GS is all the more important.

### 8. netGS step 1 on maize (drought)
*See `output/netgs_summary_DS_paper.json` and `fig_netgs_fluxes.png` (filled after the run).*

---

## Questions for Anika
1. **Objective (finding 5):** which formulation was used for the published results, the paper's Eq. 2 or the code version?
2. **Model:** OK to continue with the curated C4GEM (gap-fill + direction rules), or does the group have a curated two-cell maize model?
3. **Anchors:** OK to use Shu et al. 2022 for amino acids? Any lab data on absolute leaf metabolite levels in maize (MPI-MP)?
4. **Two-pool BOF:** is splitting bound (fixed) and free (data-driven) acceptable?
5. **WW growth data:** ask the Zhang authors (Fernie at MPI-MP) for per-line WW biomass?
6. **Expectations (finding 7):** given the weak growth signal in the BOF metabolites, should we also test adding a physiological drought constraint (CO₂ uptake / stomatal closure)?
7. **"One set of data":** what did you mean?

## Next steps (proposal)
1. Fill the remaining anchors.
2. Build the two-cell (mesophyll/bundle sheath) model on the curated network.
3. Get genotypes from maizego.org and run netGS steps 2–4 (rrBLUP on fluxes → biomass prediction) vs. classical GS, under DS.
4. Add WW if data arrive, then the G×E analysis of fluxes.

**Sources:** Zhang et al. 2021 Genome Biol 22:260 · Dal'Molin et al. 2010 Plant Physiol 154:1871 · Tong, Küken & Nikoloski 2020 Nat Commun 11:2410 (code: github.com/Hao-Tong/netGS) · Shu et al. 2022 Front Plant Sci 13:885794
