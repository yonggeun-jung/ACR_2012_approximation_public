/*
  01_build_master.do
  Build master datasets used by empirical scripts.

  Outputs
    ../data/master_gravity_agg.dta
    ../data/master_gravity_sec.dta

  Manual raw-data requirements (not shipped in replication package)
    CEPII Gravity: https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=8
      - Gravity_V202211.dta
    CEPII BACI: https://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=37
      - BACI_HS02_Y2019_V202601.csv
      - country_codes_V202601.csv
*/

clear all
set more off

* Set paths before running.
global rawdata "../data/raw"
global bacidir "../data/raw"

use "$rawdata/Gravity_V202211.dta", clear

keep year country_id_o country_id_d iso3_o iso3_d ///
     tradeflow_baci dist contig comlang_off ///
     country_exists_o country_exists_d gdp_o gdp_d

keep if year >= 2015 & year <= 2019
keep if country_exists_o == 1 & country_exists_d == 1
drop if country_id_o == country_id_d
drop if missing(dist) | dist == 0
drop if missing(tradeflow_baci)

* Select top-30 exporters by average GDP.
preserve
    collapse (mean) gdp_o, by(iso3_o)
    gsort -gdp_o
    gen rank = _n
    keep if rank <= 30
    keep iso3_o
    tempfile top
    save `top'
restore
merge m:1 iso3_o using `top', keep(match) nogen

preserve
    use `top', clear
    rename iso3_o iso3_d
    tempfile top_d
    save `top_d'
restore
merge m:1 iso3_d using `top_d', keep(match) nogen

gen ln_dist  = ln(dist)
gen ln_dist2 = ln_dist^2
egen pair_id = group(country_id_o country_id_d)

compress
save "../data/master_gravity_agg.dta", replace
di as txt "[OK] aggregate master: ../data/master_gravity_agg.dta  (N=`=_N')"

* Build sectoral bilateral panel (BACI HS02, 2019 only).
preserve
    keep if year == 2019
    keep iso3_o iso3_d dist ln_dist
    duplicates drop iso3_o iso3_d, force
    gen ln_dist2 = ln_dist^2
    tempfile dist_data
    save `dist_data'

    keep iso3_o
    rename iso3_o iso3
    duplicates drop
    tempfile top30
    save `top30'
restore

import delimited "$bacidir/country_codes_V202601.csv", clear varnames(1)
rename country_code ccode
rename country_iso3 iso3
keep ccode iso3
merge m:1 iso3 using `top30', keep(match) nogen
tempfile cw30
save `cw30'

import delimited "$bacidir/BACI_HS02_Y2019_V202601.csv", clear varnames(1)

rename i ccode
merge m:1 ccode using `cw30', keep(match) nogen
rename iso3 iso3_o
rename ccode ccode_o

rename j ccode
merge m:1 ccode using `cw30', keep(match) nogen
rename iso3 iso3_d
drop ccode

drop if iso3_o == iso3_d

capture confirm numeric variable k
if !_rc {
    tostring k, gen(k_str) format(%06.0f)
}
else {
    gen str12 k_str = k
}
replace k_str = trim(k_str)
gen hs2 = real(substr(k_str, 1, 2))
drop if missing(hs2)
drop k_str

collapse (sum) v, by(hs2 iso3_o iso3_d)
rename v trade

merge m:1 iso3_o iso3_d using `dist_data', keep(match) nogen
egen pair_id = group(iso3_o iso3_d)
compress
save "../data/master_gravity_sec.dta", replace
di as txt "[OK] sectoral master:  ../data/master_gravity_sec.dta  (N=`=_N')"
