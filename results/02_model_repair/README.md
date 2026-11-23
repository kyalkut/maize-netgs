# 2. Repairing the model

Scripts (in `code/`): `model_curation.py`

- `gapfill_reactions.json`: the 46 added reactions (43 from Dal'Molin et al. Table S1, plus 3).
- `broken_reactions.csv`: the 12 reactions switched off, with the reason for each.
- `model_curation_log.csv` (columns `reaction`, `change`): every change made, including the 28 directions set by enzyme class.
- `run_log.txt`, step 1: counts and the checks (no adenosine triphosphate or biomass without light or carbon dioxide).
