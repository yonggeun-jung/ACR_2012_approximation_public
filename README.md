# Replication Package

This repository contains code and processed replication datasets for:
**When Are Sufficient Statistics Sufficient? Bounding Gains-from-Trade Errors**.

## Repository structure

```text
data/
  master_gravity_agg.dta   Built by scripts/01_build_master.do (aggregate panel)
  master_gravity_sec.dta   Built by scripts/01_build_master.do (HS2 sector panel, 2019)
  raw/                     Large third-party raw files (git-ignored)

scripts/
  01_build_master.do       Builds master datasets from CEPII Gravity + BACI
  02_mc_validation.py      Runs Monte Carlo validation and writes figures/tables
  03_empirical.do          Runs aggregate/sectoral PPML diagnostics
  distributions.py         DGP helper functions for Monte Carlo

output/
  tables/                  CSV outputs
  figures/                 PDF outputs
```

## Requirements

- Stata 17+ (`ppmlhdfe`, `reghdfe`, `ftools`; installable from within the do-file if needed)
- Python 3.10+ with `numpy`, `scipy`, `matplotlib`

## Raw data (not distributed here)

Raw CEPII/BACI files are third-party data and are intentionally excluded from this public repository.

1. Download CEPII Gravity from [CEPII Gravity](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8) and prepare `Gravity_V202211.dta`.
2. Download CEPII BACI from [CEPII BACI](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37) and place these files in `data/raw/`:
   - `BACI_HS02_Y2019_V202601.csv`
   - `country_codes_V202601.csv`
3. In `scripts/01_build_master.do`, set:
   - `global rawdata` to the directory containing `Gravity_V202211.dta`
   - `global bacidir` to the directory containing BACI CSV files (default: `../data/raw`)

## How to run

From `scripts/`:

```bash
stata-mp -b do 01_build_master.do
python3 02_mc_validation.py
stata-mp -b do 03_empirical.do
```

If needed, use `stata` or `stata-se` instead of `stata-mp`.
