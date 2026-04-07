# Replication

**When Are Sufficient Statistics Sufficient? Bounding Gains-from-Trade Errors**  
Yonggeun Jung, University of Florida

## Layout

```
data/
  master_gravity_agg.dta   Built by scripts/01_build_master.do (CEPII aggregate panel)
  master_gravity_sec.dta   Built by scripts/01_build_master.do (BACI 2019 sector panel)
  raw/                     Large raw files (git-ignored)

scripts/
  01_build_master.do       Builds master_gravity_agg.dta + master_gravity_sec.dta
  02_mc_validation.py      Monte Carlo, CSV tables, fig_elasticity_profile_ek.pdf, fig_elasticity_profile_melitz.pdf
  03_empirical.do          Aggregate + sectoral PPML diagnostics → CSV tables
  distributions.py         DGPs for 02_mc_validation.py

output/
  tables/                  CSV
  figures/                 PDF (e.g. fig_elasticity_profile_ek.pdf, fig_elasticity_profile_melitz.pdf)
```

## Software

- Stata 17+ with `ppmlhdfe`, `reghdfe`, `ftools` (03_empirical.do can install missing packages)
- Python 3.10+ with NumPy, SciPy, Matplotlib

## Raw data setup (not in this repository)

`01_build_master.do` does **not** fetch data automatically.
You must prepare two raw sources locally before running it:

1. Download CEPII Gravity from [CEPII Gravity](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8), including **`Gravity_V202211.dta`**.
2. Download BACI from [CEPII BACI](https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37), then put **`BACI_HS02_Y2019_V202601.csv`** and **`country_codes_V202601.csv`** in `data/raw/`.
3. In `scripts/01_build_master.do`, set `global rawdata` to the folder containing `Gravity_V202211.dta`. (`global bacidir` defaults to `../data/raw`.)

The replication package ships code and generated master files only; CEPII/BACI raw files are third-party and must be downloaded manually.

## Run (from `scripts/`)

```bash
cd scripts
stata-mp -b do 01_build_master.do
python3 02_mc_validation.py
stata-mp -b do 03_empirical.do
```

Use `stata-se` or `stata` instead of `stata-mp` if needed.  
`01_build_master.do` will fail unless CEPII/BACI raw files are prepared as above.  
Steps 2–3 only need `data/master_gravity_agg.dta` and `data/master_gravity_sec.dta`.
