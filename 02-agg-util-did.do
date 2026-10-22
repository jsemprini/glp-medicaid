*==============================================================================
* 02-agg-util-did.do
*
* Aggregate state-quarter DiD (Borusyak-Jaravel-Spiess imputation, did_imputation) for the four
* GLP-1 utilization outcomes, under four control specifications, for the main and extended samples.
*
* Outcomes (state-quarter file GLP1_state_quarter_analysis.dta, via state_controls_yq.dta)
*   rx_ob_1k_adult      obesity GLP-1 prescriptions per 1,000 adult Medicaid enrollees
*   spend_ob_adult      obesity GLP-1 spending per adult enrollee
*   rx_nonob_1k_adult   non-obesity GLP-1 prescriptions per 1,000 adult enrollees   (placebo-type outcome)
*   spend_nonob_adult   non-obesity GLP-1 spending per adult enrollee                 (placebo-type outcome)
*
* Specifications (cumulative; state and quarter fixed effects in every one)
*   S1 unadjusted                       no state-level controls
*   S2 adjusted for unwinding           + ln_unwind
*   S3 + enrollment                     + L1_ln_adult  (log adult Medicaid enrollment, one-quarter lag)
*   S4 + policy controls ALL            + gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share
*                                         post_x_base_exparte3 ha_usds_assist
* With fixsample = 1 (default) every specification is estimated on the same state-quarters
* (the ones with every S4 control observed), so differences across columns come from the controls only.
* That drops DC (no NCSL politics), Nebraska (no partisan legislature) and 2021Q1 (the enrollment lag) from every
* column, so S1 here will not reproduce an unadjusted run on all 39 states and 16 quarters; set fixsample = 0 for that.
*
* Samples (the repo sample drops are kept as they were)
*   main  2021Q1-2024Q4   drops FIPS 1 12 13 20 28 37 45 46 47 48 55 56
*   ext   2021Q1-2025Q4   same drops plus FIPS 6 32 28 (no 2025 BRFSS interviews in CA and NV)
*
* For each outcome and specification: static (single tau) and dynamic (tau0-tau5, pre1-pre6,
* minn(7)) estimates, plus the joint test that pre1-pre6 are all zero.
*
* VERSION 2 (convergence patch). Version 1 used did_imputation with tol(0.001) and the default maxit(100); 23 of the 64
* fits (all of S4, and some S2/S3 dynamic fits) stopped with "Convergence of standard errors is not achieved".
* Each fit is now tried down a ladder until it returns standard errors:
*     level 1  tol(0.001) maxit(100)     identical to version 1, so every fit that converged before is reproduced exactly
*     level 2  tol(0.001) maxit(1000)
*     level 3  tol(0.01)  maxit(1000)
*     level 4  tol(0.05)  maxit(2000)
*     level 5  nose                      point estimates only (no standard errors); the table shows n/c
*     level 0  nothing worked
* The level each fit needed is saved as extra rows (term = conv_level, b = level, stat = e(Niter)), and the joint pre-trend
* statistic that did_imputation itself reports is saved as term = pretrend_native (p = e(pre_p), stat = e(pre_F)).
* Report the level for any fit above 1 in the manuscript (a looser tolerance gives less exact standard errors).
*
* Needs: did_imputation (ssc install did_imputation), state_controls_yq.dta from 01-build-state-controls.do
* Output (in $res): agg_util_did24.log / agg_util_did25.log, *_results.dta, *_long.csv, *_wide.csv
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
*----------------------------------------------------------------------------
global ctrl "$root/1-stata/controls"
global res  "$root/6-final/ii-results"
global lad1 "tol(0.001) maxit(100)"
global lad2 "tol(0.001) maxit(1000)"
global lad3 "tol(0.01) maxit(1000)"
global lad4 "tol(0.05) maxit(2000)"
global nlad 4
global enrvar L1_ln_adult
global pol_all gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist
global outcomes rx_ob_1k_adult spend_ob_adult rx_nonob_1k_adult spend_nonob_adult
global cs_S1 ""
global cs_S2 "ln_unwind"
global cs_S3 "ln_unwind $enrvar"
global cs_S4 "ln_unwind $enrvar $pol_all"
global lab_S1 "unadjusted"
global lab_S2 "adjusted for unwinding"
global lab_S3 "unwinding + enrollment"
global lab_S4 "unwinding + enrollment + policy controls ALL"

*---- helper programs (identical in 02, 03 and 04) ----------------------------
* post_est   posts one coefficient from the estimation in memory (missing if the term is absent)
* post_joint posts the joint test that the listed coefficients are all zero (uses the ones that exist)
* post_fs    posts the first-stage coefficient on glp and the cluster-robust F for it
* post_conv, post_native: see below (agg only; harmless elsewhere)
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

* post_conv   posts which rung of the convergence ladder a fit needed (term conv_level: b = level, stat = e(Niter))
* post_native posts the joint pre-trend test that did_imputation itself reports (term pretrend_native: p = e(pre_p), stat = e(pre_F))
capture program drop post_conv
program define post_conv
    args h run outc spec model lvl
    local it = .
    capture local it = e(Niter)
    post `h' ("`run'") ("`outc'") ("`spec'") ("`model'") ("conv_level") (`lvl') (.) (.) (`it') (e(N))
end

capture program drop post_native
program define post_native
    args h run outc spec model
    local pp = .
    local pf = .
    capture local pp = e(pre_p)
    capture local pf = e(pre_F)
    post `h' ("`run'") ("`outc'") ("`spec'") ("`model'") ("pretrend_native") (.) (.) (`pp') (`pf') (e(N))
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
    replace cell = "native p=" + strtrim(string(p, "%6.3f")) if term == "pretrend_native" & !missing(p)
    replace cell = "level " + strtrim(string(b, "%1.0f")) if term == "conv_level" & !missing(b)
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
        gen double d41 = bS4 - bS1 if !inlist(term, "pretrend_joint", "pretrend_native", "conv_level")
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

    use "$ctrl/state_controls_yq.dta", clear

    *---- sample ----------------------------------------------------------------
    drop if inlist(_state, 1, 12, 13, 20, 28, 37, 45, 46, 47, 48) | inlist(_state, 55, 56)
    if "`run'" == "main" {
        drop if yq_date > tq(2024q4)
    }
    if "`run'" == "ext" {
        drop if inlist(_state, 6, 32, 28)
    }
    local tag = cond("`run'" == "main", "24", "25")

    *---- treatment timing (as in the repo) ----------------------------------------
    gen time = yq
    capture drop first_tx
    bysort _state: egen first_tx = min(cond(post == 1, time, .))
    summarize time, meanonly
    local tmax = r(max)
    replace first_tx = `tmax' + 1 if missing(first_tx)
    isid _state time
    assert post == (time >= first_tx)
    bysort _state (time): assert post >= post[_n-1] if _n > 1
    xtset _state time

    *---- log and results file --------------------------------------------------------
    capture log close
    log using "$res/agg_util_did`tag'.log", replace text
    di as text _n "=== Aggregate utilization DiD, `run' sample, run on $S_DATE $S_TIME ==="
    di as text "fixsample = $fixsample   (1 = the same state-quarters in S1-S4)"
    quietly egen byte _tg = tag(_state)
    quietly count if _tg
    di as text "states in the sample before the control restriction: " r(N)
    quietly count if post == 1
    di as text "treated state-quarters (post = 1): " r(N)
    drop _tg
    foreach s of global specs {
        di as text "  `s' = ${lab_`s'}   controls: ${cs_`s'}"
    }

    tempname h
    postfile `h' str5 run str24 outcome str3 spec str12 model str24 term double b double se double p double stat long n ///
        using "$res/agg_util_did`tag'_results.dta", replace

    foreach y of global outcomes {
        foreach s of global specs {

            local cs "${cs_`s'}"
            local cm `y' `cs'
            if $fixsample == 1 {
                local cm `y' ${cs_S4}
            }
            capture drop _nm
            egen byte _nm = rowmiss(`cm')
            quietly count if _nm == 0
            local nobs = r(N)
            local ctl ""
            if "`cs'" != "" {
                local ctl "controls(`cs')"
            }

            di as text _n "-------- `run' | `y' | `s' (${lab_`s'}) | state-quarters used: `nobs' --------"

            * static: single tau
            di as text "static (single tau)"
            ereturn clear
            local lvl = 0
            forvalues l = 1/$nlad {
                if `lvl' == 0 {
                    ereturn clear
                    capture noisily did_imputation `y' _state time first_tx if _nm == 0, cluster(_state) fe(_state time) `ctl' ${lad`l'}
                    if _rc == 0 local lvl = `l'
                }
            }
            if `lvl' == 0 {
                ereturn clear
                capture noisily did_imputation `y' _state time first_tx if _nm == 0, cluster(_state) fe(_state time) `ctl' nose
                if _rc == 0 local lvl = 5
            }
            di as text "static: convergence level `lvl' (1 = as in version 1, 5 = nose, 0 = failed)"
            post_est `h' "`run'" "`y'" "`s'" "static" "tau"
            post_conv `h' "`run'" "`y'" "`s'" "static" `lvl'

            * dynamic: horizons 0-5, pretrends 1-6
            di as text "dynamic (tau0-tau5, pre1-pre6)"
            local lvl = 0
            forvalues l = 1/$nlad {
                if `lvl' == 0 {
                    ereturn clear
                    capture noisily did_imputation `y' _state time first_tx if _nm == 0, cluster(_state) fe(_state time) `ctl' pretrends(6) horizons(0/5) minn(7) ${lad`l'}
                    if _rc == 0 local lvl = `l'
                }
            }
            if `lvl' == 0 {
                ereturn clear
                capture noisily did_imputation `y' _state time first_tx if _nm == 0, cluster(_state) fe(_state time) `ctl' pretrends(6) horizons(0/5) minn(7) nose
                if _rc == 0 local lvl = 5
            }
            di as text "dynamic: convergence level `lvl' (1 = as in version 1, 5 = nose, 0 = failed)"
            forvalues k = 0/5 {
                post_est `h' "`run'" "`y'" "`s'" "dynamic" "tau`k'"
            }
            forvalues k = 1/6 {
                post_est `h' "`run'" "`y'" "`s'" "dynamic" "pre`k'"
            }
            post_joint `h' "`run'" "`y'" "`s'" "dynamic" "pre1 pre2 pre3 pre4 pre5 pre6"
            post_conv `h' "`run'" "`y'" "`s'" "dynamic" `lvl'
            post_native `h' "`run'" "`y'" "`s'" "dynamic"
        }
    }
    postclose `h'

    *---- comparison tables ---------------------------------------------------------
    di as text _n "================= `run' sample: comparison of S1-S4, aggregate utilization DiD ================="
    global nfilter `"(term == "tau" & model == "static") | (term == "tau0" & model == "dynamic")"'
    capture noisily cmp_table "$res/agg_util_did`tag'_results.dta" "$res/agg_util_did`tag'"

    log close
}
