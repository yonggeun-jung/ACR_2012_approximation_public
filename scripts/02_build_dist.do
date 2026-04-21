clear all
set more off

local pwd_now "`c(pwd)'"
if regexm("`pwd_now'", "/scripts/?$") {
    global PROJROOT ".."
}
else {
    global PROJROOT "."
}

use "$PROJROOT/data/cepii/raw/Gravity_V202211.dta", clear

keep if year == 2015
keep iso3_o iso3_d dist contig comlang_off col_dep_ever
drop if missing(dist) | dist == 0

local cp_iso "ARG AUS AUT BRA CAN CHL CHN DNK FIN FRA DEU GRC HUN IND IDN IRL ITA JPN KOR MEX NLD NZL NOR PRT ZAF ESP SWE TUR GBR USA"

gen keep_o = 0
gen keep_d = 0
foreach c of local cp_iso {
    replace keep_o = 1 if iso3_o == "`c'"
    replace keep_d = 1 if iso3_d == "`c'"
}
keep if keep_o == 1 & keep_d == 1
drop keep_o keep_d
drop if iso3_o == iso3_d

export delimited using "$PROJROOT/data/cp2015/dist_cp31.csv", replace