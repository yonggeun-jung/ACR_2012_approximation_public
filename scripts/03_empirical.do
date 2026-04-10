/*
  03_empirical.do
  Aggregate + sectoral PPML diagnostics for ACR bounds.

  Inputs (built by 01_build_master.do)
    ../data/master_gravity_agg.dta
    ../data/master_gravity_sec.dta

  Outputs
    ../output/tables/table_lambda.csv
    ../output/tables/table_empirical.csv
    ../output/tables/table_robustness_sigma.csv
    ../output/tables/table_sectoral.csv
*/

clear all
set more off

foreach pkg in reghdfe ftools ppmlhdfe {
    cap which `pkg'
    if _rc ssc install `pkg'
}
cap ftools, compile
cap reghdfe, compile

* 1) Compute domestic expenditure share lambda.

use "../data/master_gravity_agg.dta", clear
cap gen ln_dist2 = ln_dist^2

* Total exports by exporter-year
preserve
    collapse (sum) total_exports = tradeflow_baci, by(iso3_o year gdp_o)
    gen lambda = 1 - total_exports / gdp_o
    replace lambda = . if lambda <= 0 | lambda >= 1   // drop implausible
    gen ln_lambda = ln(lambda)
    
    di _n as txt "[Aggregate] Domestic expenditure share (lambda)"
    tabstat lambda, by(year) stat(mean median min max sd) format(%9.3f)
    
    * Year-to-year change in ln(lambda) per country.
    sort iso3_o year
    by iso3_o: gen dln_lambda = ln_lambda - ln_lambda[_n-1] if _n > 1
    gen abs_dln = abs(dln_lambda)
    
    di _n as txt "[Aggregate] |Delta ln lambda| (year-to-year)"
    tabstat abs_dln, stat(mean median p25 p75 max) format(%9.4f)
    
    * Representative |ln lambda_hat|: median across country-year observations.
    qui sum abs_dln, detail
    local ln_lam_med = r(p50)
    local ln_lam_p75 = r(p75)
    di as txt "Median |Delta ln lambda| = `ln_lam_med'"
    di as txt "P75    |Delta ln lambda| = `ln_lam_p75'"
    
    * Save lambda table
    collapse (mean) lambda (mean) abs_dln, by(year)
    rename lambda mean_lambda
    rename abs_dln mean_abs_dln
    export delimited using "../output/tables/table_lambda.csv", replace
restore

* Use median |Delta ln lambda| from data; fallback if missing.
cap local ln_lam = `ln_lam_med'
if missing(`ln_lam') local ln_lam = 0.05

di _n as txt "Using |ln lambda_hat| = `ln_lam'"


* 2) Aggregate gravity estimation by year.

qui sum ln_dist
local ln_d_min = r(min)
local ln_d_max = r(max)

levelsof year, local(years)
local nyrs : word count `years'
mat R = J(`nyrs', 8, .)

local row = 0
foreach y of local years {
    local ++row
    mat R[`row', 1] = `y'

    qui ppmlhdfe tradeflow_baci ln_dist if year == `y', ///
        absorb(iso3_o iso3_d) cluster(pair_id)
    mat R[`row', 2] = -_b[ln_dist]

    qui ppmlhdfe tradeflow_baci ln_dist ln_dist2 if year == `y', ///
        absorb(iso3_o iso3_d) cluster(pair_id)

    local b1 = _b[ln_dist]
    local b2 = _b[ln_dist2]
    mat R[`row', 3] = `b2'
    mat R[`row', 4] = _se[ln_dist2]
    mat R[`row', 5] = 2 * normal(-abs(`b2'/_se[ln_dist2]))

    local ea = -(`b1' + 2*`b2'*`ln_d_min')
    local eb = -(`b1' + 2*`b2'*`ln_d_max')
    mat R[`row', 6] = min(`ea', `eb')
    mat R[`row', 7] = max(`ea', `eb')
    mat R[`row', 8] = max(`ea', `eb') - min(`ea', `eb')
}

mat colnames R = year eps_pool beta2 se_b2 p_b2 eps_min eps_max Delta_eps

* 3) Bounds at sigma = {3, 4, 5}.

clear
svmat R, names(col)

foreach s in 3 4 5 {
    gen bnd_ek_s`s'  = (Delta_eps / (`s' - 1)) * `ln_lam' * 100
    gen bnd_mel_s`s' = (Delta_eps / ((`s' - 1) * eps_pool)) * `ln_lam' * 100
}

format eps_* beta2 Delta_eps %9.3f
format se_b2 p_b2 %9.4f
format bnd_* %9.2f

* 4) Display and export aggregate outputs.
di _n as txt "[Aggregate] ACR diagnostic (|ln lambda_hat| = `ln_lam')"

di _n as txt "[Aggregate] Main results (sigma = 4)"
list year eps_pool beta2 p_b2 eps_min eps_max Delta_eps ///
     bnd_ek_s4 bnd_mel_s4, noobs sep(0)

* Main table
preserve
    keep year eps_pool beta2 se_b2 p_b2 eps_min eps_max Delta_eps ///
         bnd_ek_s4 bnd_mel_s4
    rename bnd_ek_s4  bound_ek
    rename bnd_mel_s4 bound_melitz
    export delimited using "../output/tables/table_empirical.csv", replace
restore

* Sigma robustness table.
di _n as txt "[Aggregate] Sigma robustness"
list year Delta_eps bnd_mel_s3 bnd_mel_s4 bnd_mel_s5, noobs sep(0)

preserve
    keep year Delta_eps bnd_mel_s3 bnd_mel_s4 bnd_mel_s5 ///
                        bnd_ek_s3  bnd_ek_s4  bnd_ek_s5
    export delimited using "../output/tables/table_robustness_sigma.csv", replace
restore

* 5) Sectoral estimation (HS02, 2019) from prebuilt sectoral master.

use "../data/master_gravity_sec.dta", clear

qui sum ln_dist
local ln_d_min_s = r(min)
local ln_d_max_s = r(max)

levelsof hs2, local(sectors)
local nsec : word count `sectors'
mat RS = J(`nsec', 9, .)

local row = 0
foreach s of local sectors {
    local ++row
    mat RS[`row', 1] = `s'

    qui count if hs2 == `s'
    local nn = r(N)
    mat RS[`row', 2] = `nn'
    if `nn' < 100 continue

    cap qui ppmlhdfe trade ln_dist if hs2 == `s', ///
        absorb(iso3_o iso3_d) cluster(pair_id)
    if _rc continue
    mat RS[`row', 3] = -_b[ln_dist]

    cap qui ppmlhdfe trade ln_dist ln_dist2 if hs2 == `s', ///
        absorb(iso3_o iso3_d) cluster(pair_id)
    if _rc continue

    local b1s = _b[ln_dist]
    local b2s = _b[ln_dist2]
    local se2s = _se[ln_dist2]
    mat RS[`row', 4] = `b2s'
    mat RS[`row', 5] = `se2s'
    mat RS[`row', 6] = 2 * normal(-abs(`b2s'/`se2s'))

    local eas = -(`b1s' + 2*`b2s'*`ln_d_min_s')
    local ebs = -(`b1s' + 2*`b2s'*`ln_d_max_s')
    local Depss = abs(`eas' - `ebs')
    mat RS[`row', 7] = `Depss'

    mat RS[`row', 8] = (`Depss'/(4-1)) * `ln_lam' * 100

    * Melitz bound: use pooled eps_pool (linear PPML), not quadratic beta1.
    mat RS[`row', 9] = (`Depss'/((4-1)*RS[`row', 3])) * `ln_lam' * 100
}

mat colnames RS = hs2 n_obs eps_pool beta2 se_b2 p_b2 Delta_eps bnd_ek bnd_mel

clear
svmat RS, names(col)
drop if missing(eps_pool)

format eps_pool beta2 Delta_eps %9.3f
format se_b2 p_b2 %9.4f
format bnd_* %9.3f

gen sig5  = (p_b2 < 0.05)
gen sig10 = (p_b2 < 0.10)

di _n as txt "[Sectoral] Estimated sectors: `=_N'"
qui sum sig5
di as txt "[Sectoral] Significant at 5%:  `=r(sum)' / `=_N'"
qui sum sig10
di as txt "[Sectoral] Significant at 10%: `=r(sum)' / `=_N'"

gsort p_b2
list hs2 n_obs eps_pool beta2 p_b2 Delta_eps bnd_ek bnd_mel in 1/15, noobs

export delimited using "../output/tables/table_sectoral.csv", replace

di _n as txt "[Done] Tables saved to ../output/tables/"
