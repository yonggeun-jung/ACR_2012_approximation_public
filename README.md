# Replication

**When Are Sufficient Statistics Sufficient? Bounding Gains-from-Trade Formulas with an Application to NAFTA**  
Yonggeun Jung, University of Florida

## Layout

```
data/
  cepii/
    raw/                              CEPII Gravity + BACI CSVs (large; often git-ignored)
    master_gravity_agg.dta            Built by scripts/04_app_build_master_for_empiric.do
    master_gravity_sec.dta            Built by scripts/04_app_build_master_for_empiric.do
  cp2015/                             Caliendo–Parro (2015) replication inputs + dist_cp31.csv

scripts/
  01_mc_validation.py                 Monte Carlo → output/tables + output/figures
  02_build_dist.do                    Gravity → data/cp2015/dist_cp31.csv (CP country pairs)
  03_cp2015_application.py            NAFTA application → output/tables (CP bounds)
  04_app_build_master_for_empiric.do  CEPII/BACI → data/cepii/master_gravity_*.dta
  05_app_empiric.do                   PPML diagnostics → output/tables/*.csv
  acr_bound.py                        Reusable single-function ACR diagnostic
  distributions.py                    DGPs for 01_mc_validation.py

output/
  tables/                             Results tables
  figures/                            Figures from Monte Carlo
```

## Software

- Stata 17+ with `ppmlhdfe`, `reghdfe`, `ftools` (`05_app_empiric.do` can install missing packages)
- Python 3.10+ with NumPy, SciPy, Matplotlib; NAFTA application additionally uses `statsmodels`

## Raw data (not in this repository)

`04_app_build_master_for_empiric.do` does not download data. Place the following in **`data/cepii/raw/`** before running it:

1. **CEPII Gravity** — [CEPII Gravity](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8): `Gravity_V202211.dta`
2. **CEPII BACI** — [CEPII BACI](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37): `BACI_HS02_Y2019_V202601.csv`, `country_codes_V202601.csv`

The replication package ships code and (optionally) built intermediates; some third-party raw files must be obtained separately.

## Run order

From the repository root (or from `scripts/`; the Stata scripts set `PROJROOT` so paths resolve in both cases):

```bash
# 1) Build gravity masters (requires data/cepii/raw/ as above)
stata-mp -b do scripts/04_app_build_master_for_empiric.do

# 2) Empirical PPML tables (reads data/cepii/master_gravity_*.dta)
stata-mp -b do scripts/05_app_empiric.do

# 3) CP bilateral distances for Python NAFTA exercise (reads same Gravity file)
stata-mp -b do scripts/02_build_dist.do

# 4) Monte Carlo (independent of Stata masters once packages are installed)
python3 scripts/01_mc_validation.py

# 5) NAFTA / CP (2015) bounds (reads data/cp2015/* and output from step 3)
python3 scripts/03_cp2015_application.py
```

Use `stata-se` or `stata` instead of `stata-mp` if needed.

- Step 2 only needs `data/cepii/master_gravity_agg.dta` and `data/cepii/master_gravity_sec.dta`.
- Step 3 needs `data/cepii/raw/Gravity_V202211.dta` and writes `data/cp2015/dist_cp31.csv`.
- Step 5 needs CP replication inputs under `data/cp2015/` (see docstring in `scripts/03_cp2015_application.py`).

## Reusable diagnostic function

If you want to apply the bound to your own data without running the full
replication pipeline, use `scripts/acr_bound.py`:

```python
from scripts.acr_bound import compute_acr_bound

out = compute_acr_bound(
    trade_flow=trade,            # bilateral trade flow vector
    ln_distance=ln_dist,         # log distance vector
    sigma=4.0,
    ln_lambda_hat_abs=0.05,      # |ln lambda_hat|
    bilateral_controls=controls, # optional 2D array
    exporter_ids=exp_id,         # optional FE identifiers
    importer_ids=imp_id,         # optional FE identifiers
)

print(out["bound_ek"], out["bound_melitz"])
```
