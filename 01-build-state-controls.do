*==============================================================================
* 01-build-state-controls.do
*
* Merges the policy variables into the state-quarter file that every analysis file uses.
* Run once, after build_politics_panel.do, build_unwinding_panel.do and merge_state_controls.do
* (they create controls_state_yq.dta from the NCSL, CMS, KFF and Health Affairs data).
*
* Inputs
*   $prep/GLP1_state_quarter_analysis.dta   state-quarter utilization (rx_ob_1k_adult, spend_ob_adult,
*                                           rx_nonob_1k_adult, spend_nonob_adult), enroll_adult, post
*   $prep/b-data/unwinding.dta              state-quarter unwinding disenrollments ("value")
*   $ctrl/controls_state_yq.dta             governor party, % Dem senate and house, CMS unwinding measures,
*                                           Health Affairs intervention flag (policy controls)
* Output
*   $ctrl/state_controls_yq.dta   one row per state and quarter, keys: _state, yq (2021Q1 = 1) and yq_date (%tq)
*
* Controls used by the analysis files (specifications, cumulative):
*   S1 unadjusted            no state-level controls (individual models still carry the individual-level $controls)
*   S2 + unwinding           ln_unwind = log(unwinding disenrollments + 1)
*   S3 + enrollment          L1_ln_adult = log adult Medicaid enrollment, lagged one quarter
*   S4 + policy controls ALL gov_dem pct_dem_senate pct_dem_house   politics (NCSL)
*                            unw_term_pb held_proc_share           CMS unwinding intensity and procedural-termination holds
*                            post_x_base_exparte3                  early ex parte share x post-April-2023 (state capacity)
*                            ha_usds_assist                        Health Affairs intervention states (CA, NY, SC, WI)
* Enrollment is lagged because coverage can change enrollment (a bad control if contemporaneous).
*
* The log shows residualized correlations and VIFs for the full control set in the main and extended samples.
* Written for the repo's variable names (commit 6ab5f9e). Not run: Stata was not available when this was written.
*==============================================================================
clear all
set more off
set linesize 200

*---- edit these (or set once in 00-run-all.do) --------------------------------
if "$root" == "" global root "C:/Users/jsemprini/OneDrive - Des Moines University/9-Internal-Collab/a-Devine_GLP1"
global prep "$root/5-iv-prep"
global ctrl "$root/1-stata/controls"
global res  "$root/6-final/ii-results"
global pol_all gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist
*----------------------------------------------------------------------------

* every policy variable brought in by the merge
global newvars gov_dem pct_dem_senate pct_dem_house pct_dem2_senate pct_dem2_house nonpartisan_leg ///
               unw_term_pb unw_proc_pb unw_cum_term_pb unw_rate_term held_proc_share base_exparte3 ///
               post_x_base_exparte3 cum_disenroll cum_exparte_renewed cum_proc_share ha_usds_assist ///
               cms_partial no_politics

capture program drop vif_resid
program define vif_resid
    * SD, residualized SD and VIF of each variable after removing state and quarter effects
    syntax varlist
    foreach v of local varlist {
        capture drop rr_`v'
        quietly areg `v' i.yq, absorb(_state)
        quietly predict double rr_`v', residuals
    }
    display as text _n %-26s "variable" %12s "SD raw" %14s "SD after FE" %10s "VIF"
    foreach v of local varlist {
        local others : list varlist - v
        local rhs ""
        foreach o of local others {
            local rhs `rhs' rr_`o'
        }
        quietly summarize `v'
        local sd0 = r(sd)
        quietly summarize rr_`v'
        local sd1 = r(sd)
        local vif = .
        capture quietly regress rr_`v' `rhs'
        if _rc == 0 local vif = 1/(1 - e(r2))
        local flag = cond(missing(`vif'), "  (not computed)", cond(`vif' > 10, "  <<< VIF above 10", cond(`vif' > 5, "  < VIF above 5", "")))
        display as text %-26s "`v'" %12.4f `sd0' %14.4f `sd1' %10.2f `vif' "`flag'"
    }
    display as text _n "Residualized correlations (state and quarter effects removed):"
    local rlist ""
    foreach v of local varlist {
        local rlist `rlist' rr_`v'
    }
    pwcorr `rlist', obs
    drop `rlist'
end

cd "$res"
capture log close
log using "build_state_controls.log", replace text

*---- 1. state-quarter file ----------------------------------------------------
use "$prep/GLP1_state_quarter_analysis.dta", clear
capture confirm variable yq
if _rc == 0 rename yq full_yq
foreach v in year quarter {
    capture confirm numeric variable `v'
    if _rc destring `v', replace force
}
egen yq = group(year quarter)
gen int yq_date = yq(year, quarter)
format yq_date %tq
* yq in unwinding.dta and in the individual-level merges is egen group(year quarter) with 2021Q1 = 1
assert yq == yq_date - tq(2021q1) + 1
isid _state yq
foreach v in post enroll_adult rx_ob_1k_adult spend_ob_adult rx_nonob_1k_adult spend_nonob_adult {
    confirm variable `v'
}

*---- 2. unwinding and enrollment (same construction as the repo) ---------------
merge 1:1 _state yq using "$prep/b-data/unwinding.dta", generate(_mu)
tab _mu
drop if _mu == 2
drop _mu
confirm variable value
replace value = 0 if missing(value)
gen ln_unwind = log(value + 1)
gen ln_adult  = log(enroll_adult)
xtset _state yq
gen L1_ln_adult = L.ln_adult
label var ln_unwind   "log(unwinding disenrollments + 1)"
label var ln_adult    "log adult Medicaid enrollment"
label var L1_ln_adult "log adult Medicaid enrollment, lagged one quarter"

*---- 3. policy variables ----------------------------------------------------------
foreach v of global newvars {
    confirm new variable `v'
}
merge 1:1 _state yq_date using "$ctrl/controls_state_yq.dta", keepusing($newvars) keep(master match) generate(_mc)
tab _mc
* 3 = matched. 1 = a state-quarter outside the controls panel (2021Q1-2025Q4, 50 states + DC); its policy variables are missing.
drop _mc

count
sort _state yq
save "$ctrl/state_controls_yq.dta", replace

*---- 4. collinearity of the full control set in the estimation samples -------------
* main = 2021Q1-2024Q4 ; ext = through 2025Q4 without CA, NV, MS (no 2025 BRFSS interviews in CA and NV)
foreach run in main ext {
    preserve
    drop if inlist(_state, 1, 12, 13, 20, 28, 37, 45, 46, 47, 48) | inlist(_state, 55, 56)
    if "`run'" == "main" drop if yq_date > tq(2024q4)
    if "`run'" == "ext" drop if inlist(_state, 6, 32, 28)
    egen byte _nm = rowmiss(ln_unwind L1_ln_adult $pol_all)
    di _n as text "======== `run' sample: state-quarters with every control observed ========"
    tab _nm
    keep if _nm == 0
    egen byte _tag = tag(_state)
    quietly count if _tag
    di as text "states in the complete-case sample: " r(N)
    quietly count
    di as text "state-quarter cells: " r(N)
    di as text "States lost to missing controls: DC (no NCSL data) and Nebraska (no partisan legislature)."
    vif_resid ln_unwind L1_ln_adult $pol_all
    restore
}
di _n as text "Reading the table: a VIF above 10 means that control is almost a linear combination of the others once"
di as text "state and quarter effects are removed. If unw_term_pb overlaps with ln_unwind, drop one of them from pol_all."

log close
