# 01_model_curation

Repair C4GEM: +46 gap-fill reactions (43 from Dal'Molin 2010 Table S1), 12 broken reactions blocked, 28 directions set by EC rules.

**Scripts:** model_curation.py (called by bof_pipeline.py, step 1)

**Inputs:**
- `inputs/Maize_C4GEM_vs1.0.xml` (from `input/Maize_C4GEM_vs1.0.xml`)
- `inputs/gapfill_reactions.json` (from `curation/gapfill_reactions.json`)
- `inputs/broken_reactions.csv` (from `curation/broken_reactions.csv`)

**Outputs:**
- `outputs/model_curation_log.csv` (from `output/model_curation_log.csv`)
