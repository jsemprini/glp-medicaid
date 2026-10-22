*==============================================================================
* 03-ind-did.do
*
* Individual-level DiD (Borusyak-Jaravel-Spiess imputation, did_imputation) on the BRFSS microdata,
* four control specifications, main and extended samples.
*
* Outcomes
*   bmi          BMI (kg/m2) = _bmi5 / 100, among Medicaid enrollees only (medicaid == 1)
*   medicaid100  100 x Medicaid as primary coverage, whole sample   (effects are percentage points)
*
* Treatment   glp = coverage of GLP-1 anti-obesity drugs in the state in the interview quarter (as built in
*             0ii-clean-final.do). time = quarter index within the estimation sample; first_tx = first treated
*             quarter of the state (missing = never treated).
*
* Specifications (cumulative). The individual-level controls ($controls) are fixed effects in every one,
* as in 3-individual24.do, together with state and quarter fixed effects.
*   S1 unadjusted                       no state-level controls
*   S2 adjusted for unwinding           + ln_unwind
*   S3 + enrollment                     + L1_ln_adult  (log adult Medicaid enrollment, one-quarter lag)
*   S4 + policy controls ALL            + gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share
*                                         post_x_base_exparte3 ha_usds_assist
* With fixsample = 1 (default) S1-S4 use the same individuals (every S4 control observed), so the columns
* differ only by the controls. That drops DC (no NCSL politics), Nebraska (no partisan legislature) and the
* first quarter (enrollment lag) from every column.
*
* Samples
*   main  2021Q1-2024Q4   drops diabetes, FIPS 1 12 13 20 28 37 45 46 47 48 55 56
*   ext   2021Q1-2025Q4   same drops plus FIPS 6 32 28.  (3-individual25.do dropped survey_year 2025 by mistake;
*                         that line is not repeated here.)
*
* Estimates per outcome and specification: static tau; dynamic tau0-tau5 and pre1-pre6 with the joint pretrend test.
* Weights: [aw=_llcpwt]. Standard errors clustered by state.
*
* Runtime: the medicaid100 calls on the whole sample are slow (about 25 minutes per call in your earlier logs).
* Use $runs, $outsel and $specs (see 00-run-all.do) to split the work.
*
* Needs: did_imputation, state_controls_yq.dta from 01-build-state-controls.do
* Output (in $res): ind_did24.log / ind_did25.log, *_results.dta, *_long.csv, *_wide.csv
* Not run: Stata was not available when this was written.
*==============================================================================
clear all
set more off
set linesize 200

*---- edit these (or set them once in 00-run-all.do) -------------------------
if "$root"      == "" global root "C:/Users/jsemprini/OneDrive - Des Moines University/9-Internal-Collab/a-Devine_GLP1"
if "$runs"      == "" global runs main ext
if "$specs"     == "" global specs S1 S2 S3 S4
if "$fixsample" == "" global fixsample 1
if "$outsel"    == "" global outsel bmi medicaid
if "$dyn"       == "" global dyn 1
*----------------------------------------------------------------------------
global ctrl   "$root/1-stata/controls"
global res    "$root/6-final/ii-results"
global brfss  "$root/1-stata/brfss_bmi2021-2025analytic.dta"
global tolopt "tol(0.01) maxit(1000)"
global enrvar L1_ln_adult
global pol_all gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist
global controls dum_age2 dum_age3 dum_age4 dum_age5 dum_age6 dum_age7 dum_age8 dum_age9 dum_race1 dum_race2 dum_race5 dum_sex2 dum_educ2 dum_educ3 dum_educ4 married metro parent
global cs_S1 ""
global cs_S2 "ln_unwind"
global cs_S3 "ln_unwind $enrvar"
global cs_S4 "ln_unwind $enrvar $pol_all"
global lab_S1 "unadjusted"
global lab_S2 "adjusted for unwinding"
global lab_S3 "unwinding + enrollment"
global lab_S4 "unwinding + enrollment + policy controls ALL"
* outcome definitions: name, condition for the estimation sample, label in the results file
global cond_bmi         "medicaid == 1"
global cond_medicaid ""
global name_bmi         "bmi_if_medicaid"
global name_medicaid "medicaid"

*---- helper programs (identical in 02, 03 and 04) ----------------------------
* post_est   posts one coefficient from the estimation in memory (missing if the term is absent)
* post_joint posts the joint test that the listed coefficients are all zero (uses the ones that exist)
* post_fs    posts the first-stage coefficient on glp and the cluster-robust F for it
* cmp_table  builds the side-by-side comparison of the four specifications from the results file
capture program drop post_est
program define post_est
    args h run outc spec model term
    local b  = .
    local se = .
    local p  = .
    capture local b  = _b[`term']
    capture local se = _se[`term']
    if `se' < . {
        if `se' > 0 {
            capture quietly test `term'
            if _rc == 0 local p = r(p)
        }
    }
    post `h' ("`run'") ("`outc'") ("`spec'") ("`model'") ("`term'") (`b') (`se') (`p') (.) (e(N))
end

capture program drop post_joint
program define post_joint
    args h run outc spec model terms
    local p    = .
    local stat = .
    local have ""
    foreach t of local terms {
        capture local chk = _b[`t']
        if _rc == 0 local have `have' `t'
    }
    if "`have'" != "" {
        capture quietly test `have'
        if _rc == 0 {
            local p    = r(p)
            local stat = cond(missing(r(chi2)), r(F), r(chi2))
        }
    }
    post `h' ("`run'") ("`outc'") ("`spec'") ("`model'") ("pretrend_joint") (.) (.) (`p') (`stat') (e(N))
end

capture program drop post_fs
program define post_fs
    args h run outc spec x
    local b  = .
    local se = .
    local p  = .
    local F  = .
    capture local b  = _b[glp]
    capture local se = _se[glp]
    capture quietly test glp
    if _rc == 0 {
        local p = r(p)
        local F = r(F)
    }
    post `h' ("`run'") ("`outc'") ("`spec'") ("FS") ("`x'") (`b') (`se') (`p') (`F') (e(N))
end

capture program drop cmp_table
program define cmp_table
    * usage: global nfilter `"<condition selecting the rows whose N is shown>"'
    *        cmp_table "<results.dta>" "<output stem>"
    args resfile stem
    preserve
    use "`resfile'", clear
    gen long seq = _n
    capture noisily export delimited using "`stem'_long.csv", replace
    gen str60 cell = ""
    replace cell = strtrim(string(b, "%12.4f")) + cond(p < .01, "***", cond(p < .05, "**", cond(p < .1, "*", ""))) + " (" + strtrim(string(se, "%12.4f")) + ")" if !missing(b) & !missing(se)
    replace cell = cell + " F=" + strtrim(string(stat, "%8.1f")) if model == "FS" & !missing(stat) & cell != ""
    replace cell = "joint p=" + strtrim(string(p, "%6.3f")) if term == "pretrend_joint" & !missing(p)
    replace cell = "n/a" if cell == ""
    gen str14 nn = string(n, "%12.0fc")
    bysort run outcome model term: egen long ord = min(seq)
    keep run outcome model term spec cell nn b ord
    reshape wide cell nn b, i(run outcome model term) j(spec) string
    sort ord
    local cl ""
    local nl ""
    foreach s in S1 S2 S3 S4 {
        capture confirm variable cell`s'
        if _rc == 0 {
            local cl `cl' cell`s'
            local nl `nl' nn`s'
        }
    }
    local dl ""
    capture confirm variable bS1
    local ok1 = _rc
    capture confirm variable bS4
    if `ok1' == 0 & _rc == 0 {
        gen double d41 = bS4 - bS1 if term != "pretrend_joint"
        format d41 %9.4f
        local dl d41
    }
    format `cl' %-34s
    local oldls = c(linesize)
    set linesize 250
    di _n as text "Cells are coefficient (std. err.) with * p<0.10, ** p<0.05, *** p<0.01; d41 = S4 minus S1 point estimate."
    list run outcome model term `cl' `dl', noobs sepby(outcome)
    di _n as text "Observations used in each specification:"
    list run outcome model term `nl' if $nfilter, noobs sepby(outcome)
    set linesize `oldls'
    drop ord
    capture noisily export delimited using "`stem'_wide.csv", replace
    restore
end
*------------------------------------------------------------------------------

foreach run of global runs {

    use "$brfss", clear

    *---- quarter keys -------------------------------------------------------------
    * yq in the analytic file is the actual interview year-quarter as a %tq date
    rename yq yq_date
    assert inrange(yq_date, tq(2021q1), tq(2025q4))
    gen int yq = yq_date - tq(2021q1) + 1

    *---- sample (same drops as 3-individual24.do) -----------------------------------
    drop if diabetes == 1
    drop if inlist(_state, 1, 12, 13, 20, 28, 37, 45, 46, 47, 48) | inlist(_state, 55, 56)
    if "`run'" == "main" {
        drop if yq_date > tq(2024q4)
    }
    if "`run'" == "ext" {
        drop if inlist(_state, 6, 32, 28)
    }
    local tag = cond("`run'" == "main", "24", "25")

    *---- outcomes -------------------------------------------------------------------
    * _bmi5 is BMI x 100 in the BRFSS file (0ii-clean-final.do cuts at 1850, 2500, ...). Rescale and check.
    quietly summarize _bmi5
    assert r(mean) > 1500 & r(mean) < 5000
    gen double bmi = _bmi5 / 100
    label var bmi "BMI (kg/m2) = _bmi5 / 100"

    *---- state-level controls --------------------------------------------------------
    foreach v in ln_unwind $enrvar $pol_all {
        capture confirm variable `v'
        if _rc == 0 {
            di as error "`v' already exists in the BRFSS file; rename it before merging"
            exit 110
        }
    }
    merge m:1 _state yq_date using "$ctrl/state_controls_yq.dta", keepusing(ln_unwind $enrvar $pol_all) keep(master match) generate(_ms)
    tab _ms
    drop _ms

    *---- treatment timing ---------------------------------------------------------------
    egen time = group(yq_date)
    bysort _state (time): egen first_tx = min(cond(glp == 1, time, .))
    assert glp == 0 if missing(first_tx)
    assert glp == (time >= first_tx) if !missing(first_tx)

    *---- log -------------------------------------------------------------------------------
    capture log close
    log using "$res/ind_did`tag'.log", replace text
    di as text _n "=== Individual-level DiD, `run' sample, run on $S_DATE $S_TIME ==="
    di as text "fixsample = $fixsample   (1 = the same individuals in S1-S4)"
    quietly egen byte _tg = tag(_state)
    quietly count if _tg
    di as text "states in the sample: " r(N)
    drop _tg
    quietly count
    di as text "individuals in the sample (before outcome and control restrictions): " r(N)
    egen byte _nmall = rowmiss(ln_unwind $enrvar $pol_all)
    di as text "Individuals with a missing state-level control after 2021Q1, by state (2021Q1 loses the enrollment lag everywhere):"
    capture noisily tab _state if _nmall > 0 & yq_date > tq(2021q1)
    drop _nmall
    foreach s of global specs {
        di as text "  `s' = ${lab_`s'}   state-level controls: ${cs_`s'}"
    }

    tempname h
    postfile `h' str5 run str24 outcome str3 spec str12 model str24 term double b double se double p double stat long n ///
        using "$res/ind_did`tag'_results.dta", replace

    foreach y in $outsel {
        local cond "${cond_`y'}"
        local oname "${name_`y'}"
        foreach s of global specs {

            local cs "${cs_`s'}"
            local cm `y' $controls _llcpwt `cs'
            if $fixsample == 1 {
                local cm `y' $controls _llcpwt ${cs_S4}
            }
            capture drop _nm
            egen byte _nm = rowmiss(`cm')
            local ifc "_nm == 0"
            if "`cond'" != "" {
                local ifc "_nm == 0 & `cond'"
            }
            quietly count if `ifc'
            local nobs = r(N)
            local ctl ""
            if "`cs'" != "" {
                local ctl "controls(`cs')"
            }

            di as text _n "-------- `run' | `oname' | `s' (${lab_`s'}) | individuals used: `nobs' --------"

            di as text "static (single tau)"
            ereturn clear
            capture noisily did_imputation `y' _state time first_tx if `ifc' [aw = _llcpwt], cluster(_state) fe(_state time $controls) `ctl' $tolopt
            post_est `h' "`run'" "`oname'" "`s'" "static" "tau"

            if $dyn == 1 {
                di as text "dynamic (tau0-tau5, pre1-pre6)"
                ereturn clear
                capture noisily did_imputation `y' _state time first_tx if `ifc' [aw = _llcpwt], cluster(_state) fe(_state time $controls) `ctl' pretrends(6) horizons(0/5) $tolopt
                forvalues k = 0/5 {
                    post_est `h' "`run'" "`oname'" "`s'" "dynamic" "tau`k'"
                }
                forvalues k = 1/6 {
                    post_est `h' "`run'" "`oname'" "`s'" "dynamic" "pre`k'"
                }
                post_joint `h' "`run'" "`oname'" "`s'" "dynamic" "pre1 pre2 pre3 pre4 pre5 pre6"
            }
        }
    }
    postclose `h'

    *---- comparison tables -------------------------------------------------------------------
    di as text _n "================= `run' sample: comparison of S1-S4, individual DiD ================="
    global nfilter `"(term == "tau" & model == "static") | (term == "tau0" & model == "dynamic")"'
    capture noisily cmp_table "$res/ind_did`tag'_results.dta" "$res/ind_did`tag'"

    log close
}
