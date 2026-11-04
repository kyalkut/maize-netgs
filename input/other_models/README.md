# Other maize models (downloaded 4 Nov 2026, not in git)

| File | Model | Source |
|---|---|---|
| `plphys_v166_3_1659_s1.zip` → `simons2014/sbml/TS4_Model.xml` (Suppl. Table S4), `simons2014/245787supptbl3.xls` (Excel model), `simons2014/245787supptbl1.xls` (biomass measurements) | Simons et al. 2014, Plant Physiol 166:1659, two-cell leaf model, 8,525 reactions | https://doi.org/10.1104/pp.114.245787, supplementary data zip |
| `pone.0151722.s019.xml` | Bogart & Myers 2016, iEB5204 (single cell, comprehensive) | https://doi.org/10.1371/journal.pone.0151722, S1 Model |
| `pone.0151722.s021.xml` | Bogart & Myers 2016, iEB2140x2 (mesophyll + bundle sheath) | same, S3 Model |

To rebuild `simons2014/`: unzip the zip, then unzip `245787supptbl4.zip` into `simons2014/sbml/`.
Used by `compare_models.py`.
