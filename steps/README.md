# Pipeline steps

Copies of each step's inputs and outputs. Made by `make_steps.py`; the scripts still use input/, curation/ and output/.

- `01_model_curation/`: Repair C4GEM: +46 gap-fill reactions (43 from Dal'Molin 2010 Table S1), 12 broken reactions blocked, 28 directions set by EC rules.
- `02_metabolite_selection/`: 351 annotated Zhang features -> 52 name matches (50 metabolites) -> 35 producible -> 31 with clean data.
- `03_rescaling/`: log2 data -> mmol/gDW: relative factor 2^(x - WW mean), x anchor, + Table III part, x mass factor f to 1 g DW (BioModTool; Schulz 2021).
- `04_bof_feasibility_and_sensitivity/`: Max biomass of all 770 BOFs (FBA); placeholder anchors x0.1 and x10.
- `05_netgs_step1/`: Per DS line: s_model, growth rescaling from dry mass, B73 reference (pFBA), QP closest to reference.
- `06_data_analyses/`: G x E consistency per BOF metabolite (WW vs DS line correlation); ridge CV prediction of DS dry mass.
- `08_fluxsum/`: Flux-sum phi = 1/2 sum|S_ij v_j| of each measured metabolite per DS line (from step 05 fluxes); Spearman r with the measured DS log2 level across lines.
- `10_variant_without_24/`: Variant: the 24 metabolites without a literature anchor left out (BOF = Table III + 7 amino acids, 37 components, 1 g DW); BOFs, netGS step 1 (paper/code x serine pathway/RuBisCO bounds), FluxSum, data analyses; compared with the main run (61 components).
- `09_model_comparison/`: C4GEM vs Simons et al. 2014 vs Bogart & Myers 2016 (iEB5204, iEB2140x2): size, cells, genes, biomass, overlap with Zhang metabolites.
- `11_netgs_steps2-4/`: Genotypes (ZEAMAP AMP SNPs, B73 v4) -> kernel K = ZZ'; netGS step 2 (rrBLUP per flux, 50 x 3-fold CV), step 3 (QP towards predicted fluxes with the test line's BOF), step 4 (r with DS dry mass) vs classical GS on dry mass; both BOF variants.
- `07_validation_arabidopsis/`: Python port vs published netGS on the Arabidopsis example (s_model identical, 0.0035036); paper (Eq. 2) vs code objective: same flux patterns, r = 1.00 (cvxpy).
