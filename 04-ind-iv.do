*==============================================================================
* 04-ind-iv.do
*
* Individual-level IV on the BRFSS microdata. Coverage is the instrument for GLP-1 utilization and spending:
*
*   first stage   x   = a + b*glp + state FE + quarter FE + individual controls + state-level controls
*   second stage  y   = c + d*x_hat + (same controls), instrument glp, 2SLS
*
* Endogenous variable x, one regression at a time (state-quarter measures merged onto the individuals)
*   rx_ob_1k_adult      obesity GLP-1 prescriptions per 1,000 adult Medicaid enrollees
*   spend_ob_adult      obesity GLP-1 spending per adult enrollee
*   rx_nonob_1k_adult   PLACEBO: non-obesity GLP-1 prescriptions per 1,000 adult enrollees
*   spend_nonob_adult   PLACEBO: non-obesity GLP-1 spending per adult enrollee
* The coverage mandate covers anti-obesity use only, so glp should move the obesity measures (strong first stage)
* and should not move the non-obesity ones (weak or zero first stage, no usable second stage).
*
* Outcomes   bmi (BMI = _bmi5 / 100) among Medicaid enrollees (medicaid == 1);  medicaid, whole sample.
*
* This reverses the orientation used in 3-individual24iv.do, where glp was the instrumented variable and utilization
* the instrument. Here glp (coverage) is the instrument, as requested.
*
* Specifications S1-S4 are the same cumulative control sets as in 02 and 03; state and quarter fixed effects and the
* individual controls ($controls) are in every one:
*   S1 unadjusted | S2 + ln_unwind | S3 + L1_ln_adult | S4 + policy controls ALL (gov_dem pct_dem_senate pct_dem_house
*   unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist).  fixsample = 1 keeps one common sample across S1-S4.
* fixsample = 1 drops DC, Nebraska and 2021Q1 (enrollment lag) from every column; set fixsample = 0 to let each
* specification use all of its own observations.
*
* Rows in the comparison tables:  model FS = first stage (coefficient on glp, with the cluster-robust F for glp);
*                                 model 2SLS = coefficient on the instrumented x.
*
* Estimation   ivregress 2sls ... (x = glp) [pweight=_llcpwt], absorb(_state yq) vce(cluster _state)   (Stata 19)
*              first stage re-run with areg so the F statistic is the cluster-robust F for the excluded instrument.
*
* Samples   main 2021Q1-2024Q4 ;  ext 2021Q1-2025Q4 without FIPS 6 32 28 (same drops as the repo).
*
* Needs: state_controls_yq.dta from 01-build-state-controls.do
* Output (in $res): ind_iv24.log / ind_iv25.log, *_results.dta, *_long.csv, *_wide.csv
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
*----------------------------------------------------------------------------
global ctrl   "$root/1-stata/controls"
global res    "$root/6-final/ii-results"
global brfss  "$root/1-stata/brfss_bmi2021-2025analytic.dta"
global enrvar L1_ln_adult
global pol_all gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist
global controls dum_age2 dum_age3 dum_age4 dum_age5 dum_age6 dum_age7 dum_age8 dum_age9 dum_race1 dum_race2 dum_race5 dum_sex2 dum_educ2 dum_educ3 dum_educ4 married metro parent
global endog rx_ob_1k_adult spend_ob_adult rx_nonob_1k_adult spend_nonob_adult
global cs_S1 ""
global cs_S2 "ln_unwind"
global cs_S3 "ln_unwind $enrvar"
global cs_S4 "ln_unwind $enrvar $pol_all"
global lab_S1 "unadjusted"
global lab_S2 "adjusted for unwinding"
global lab_S3 "unwinding + enrollment"
global lab_S4 "unwinding + enrollment + policy controls ALL"
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
    rename yq yq_date
    assert inrange(yq_date, tq(2021q1), tq(2025q4))
    gen int yq = yq_date - tq(2021q1) + 1

    *---- sample ----------------------------------------------------------------------
    drop if diabetes == 1
    drop if inlist(_state, 1, 12, 13, 20, 28, 37, 45, 46, 47, 48) | inlist(_state, 55, 56)
    if "`run'" == "main" {
        drop if yq_date > tq(2024q4)
    }
    if "`run'" == "ext" {
        drop if inlist(_state, 6, 32, 28)
    }
    local tag = cond("`run'" == "main", "24", "25")

    *---- outcomes ------------------------------------------------------------------------
    quietly summarize _bmi5
    assert r(mean) > 1500 & r(mean) < 5000
    gen double bmi = _bmi5 / 100
    label var bmi "BMI (kg/m2) = _bmi5 / 100"
    assert inlist(glp, 0, 1)

    *---- state-quarter utilization and controls ----------------------------------------------
    foreach v in $endog ln_unwind $enrvar $pol_all {
        capture confirm variable `v'
        if _rc == 0 {
            di as error "`v' already exists in the BRFSS file; rename it before merging"
            exit 110
        }
    }
    merge m:1 _state yq_date using "$ctrl/state_controls_yq.dta", keepusing($endog ln_unwind $enrvar $pol_all) keep(master match) generate(_ms)
    tab _ms
    drop _ms

    *---- log ---------------------------------------------------------------------------------------
    capture log close
    log using "$res/ind_iv`tag'.log", replace text
    di as text _n "=== Individual-level IV, `run' sample, run on $S_DATE $S_TIME ==="
    di as text "instrument: glp (coverage).  fixsample = $fixsample   (1 = the same individuals in S1-S4)"
    quietly egen byte _tg = tag(_state)
    quietly count if _tg
    di as text "states in the sample: " r(N)
    drop _tg
    quietly count
    di as text "individuals in the sample (before outcome and control restrictions): " r(N)
    egen byte _nmall = rowmiss($endog ln_unwind $enrvar $pol_all)
    di as text "Individuals with a missing state-quarter variable after 2021Q1, by state (2021Q1 loses the enrollment lag everywhere):"
    capture noisily tab _state if _nmall > 0 & yq_date > tq(2021q1)
    drop _nmall
    foreach s of global specs {
        di as text "  `s' = ${lab_`s'}   state-level controls: ${cs_`s'}"
    }

    tempname h
    postfile `h' str5 run str24 outcome str3 spec str12 model str24 term double b double se double p double stat long n ///
        using "$res/ind_iv`tag'_results.dta", replace

    foreach y in $outsel {
        local cond "${cond_`y'}"
        local oname "${name_`y'}"
        foreach x of global endog {
            foreach s of global specs {

                local cs "${cs_`s'}"
                local cm `y' `x' glp $controls _llcpwt `cs'
                if $fixsample == 1 {
                    local cm `y' `x' glp $controls _llcpwt ${cs_S4}
                }
                capture drop _nm
                egen byte _nm = rowmiss(`cm')
                local ifc "_nm == 0"
                if "`cond'" != "" {
                    local ifc "_nm == 0 & `cond'"
                }
                quietly count if `ifc'
                local nobs = r(N)

                di as text _n "-------- `run' | `oname' | x = `x' | `s' (${lab_`s'}) | individuals used: `nobs' --------"

                * first stage: x on glp, same controls and fixed effects
                ereturn clear
                capture quietly areg `x' glp i.yq $controls `cs' [pweight = _llcpwt] if `ifc', absorb(_state) vce(cluster _state)
                local fsrc = _rc
                post_fs `h' "`run'" "`oname'" "`s'" "`x'"
                if `fsrc' == 0 {
                    quietly test glp
                    di as text "first stage:  coef on glp = " %10.4f _b[glp] "   se = " %9.4f _se[glp] "   cluster-robust F for glp = " %9.2f r(F)
                }
                if `fsrc' != 0 {
                    di as error "first stage failed (rc = `fsrc')"
                }

                * second stage: 2SLS, glp instruments x
                ereturn clear
                capture quietly ivregress 2sls `y' i.($controls) `cs' (`x' = glp) [pweight = _llcpwt] if `ifc', absorb(_state yq) vce(cluster _state)
                local ivrc = _rc
                if `ivrc' == 0 {
                    di as text "second stage: coef on `x' = " %10.4f _b[`x'] "   se = " %9.4f _se[`x']
                    capture noisily ivregress
                }
                if `ivrc' != 0 {
                    di as error "ivregress failed (rc = `ivrc')"
                }
                post_est `h' "`run'" "`oname'" "`s'" "2SLS" "`x'"
            }
        }
    }
    postclose `h'

    *---- comparison tables -----------------------------------------------------------------------------
    di as text _n "================= `run' sample: comparison of S1-S4, individual IV (glp instruments x) ================="
    global nfilter `"model == "2SLS""'
    capture noisily cmp_table "$res/ind_iv`tag'_results.dta" "$res/ind_iv`tag'"

    log close
}
