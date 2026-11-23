# Validation on Arabidopsis (sections 5 and 6)

Scripts (in `code/`): `validation/validate_arabidopsis.py`, `validation/compare_objectives.py`

Reproduces s_model = 0.0035036 of Tong et al. (2020) on their Arabidopsis data and compares Equation 2 with
the code objective (median per-reaction correlation 1.00). Needs the authors' repository:
`git clone https://github.com/Hao-Tong/netGS`, then `NETGS_REPO=netGS/netGS/ python validation/validate_arabidopsis.py` from `code/`.
