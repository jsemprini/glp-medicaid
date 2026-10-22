ESTIMATION PACKAGE: state controls merged into aggregate DiD, individual DiD and individual IV
(for the repo github.com/jsemprini/glp-medicaid at commit 6ab5f9e; Stata 19; NOT run in Stata: it was not available when this was written)

SETUP
  1. Put the state controls package files in  $root/1-stata/controls  (build_politics_panel.do, build_unwinding_panel.do,
     merge_state_controls.do and their CSV inputs) and run those three do files. They create controls_state_yq.dta.
  2. Put the files in this folder in  $root/1-stata/est  (or change $code in 00-run-all.do).
  3. ssc install did_imputation, replace ; ssc install reghdfe, replace ; ssc install ftools, replace
  4. Edit the paths and switches at the top of 00-run-all.do, then run it, or run 01-04 one at a time.

FILES
  00-run-all.do               paths, switches, calls 01-04
  01-build-state-controls.do  merges unwinding, enrollment (lagged) and the policy controls onto GLP1_state_quarter_analysis.dta;
                              writes state_controls_yq.dta; log shows residualized correlations and VIFs of the full control set
  02-agg-util-did.do          aggregate DiD, four utilization outcomes: rx_ob_1k_adult spend_ob_adult rx_nonob_1k_adult spend_nonob_adult
  03-ind-did.do               individual imputation DiD: bmi | medicaid == 1, and medicaid100 (whole sample)
  04-ind-iv.do                individual IV: glp instruments each of the four utilization and spending measures (obesity and placebo
                              non-obesity); outcomes bmi | medicaid == 1 and medicaid100

SPECIFICATIONS (cumulative, in every log)
  S1 unadjusted                       state and quarter fixed effects (individual models also carry the individual controls)
  S2 adjusted for unwinding           + ln_unwind
  S3 + enrollment                     + L1_ln_adult
  S4 + policy controls ALL            + gov_dem pct_dem_senate pct_dem_house unw_term_pb held_proc_share post_x_base_exparte3 ha_usds_assist
  fixsample = 1 puts S1-S4 on the same observations.

OUTPUT (all in $res = $root/6-final/ii-results)
  agg_util_did24.log   agg_util_did25.log     main (2021Q1-2024Q4) and extended (through 2025Q4) aggregate DiD
  ind_did24.log        ind_did25.log          individual DiD
  ind_iv24.log         ind_iv25.log           individual IV (rows FS = first stage on glp with F; 2SLS = instrumented coefficient)
  *_results.dta, *_long.csv, *_wide.csv       every estimate and the S1-S4 comparison table in machine-readable form
  Each log ends with the comparison table: coefficient (std. err.) for S1-S4 side by side, d41 = S4 minus S1,
  the joint pretrend p-value for the dynamic models, and the number of observations in each specification.

SWITCHES
  runs      main ext         which samples
  specs     S1 S2 S3 S4      which specifications
  fixsample 1                same observations in every specification
  outsel    bmi medicaid100  which individual outcomes (03 and 04)
  dyn       1                03 only: set to 0 to skip the horizons / pretrends models

RUNTIME
  The medicaid100 DiD on the whole sample took about 25 minutes per call in your earlier logs. 03 makes 8 such calls per sample
  (2 per specification), so plan for hours, or split the work with runs, specs and outsel.

KNOWN LIMITS
  - ln_unwind and unw_term_pb may overlap. Check the VIFs in build_state_controls.log; if one is above 10, drop one from pol_all.
  - DC (no NCSL politics) and Nebraska (no partisan legislature) drop out of every column when fixsample = 1.
  - The lag in L1_ln_adult drops 2021Q1 from every column when fixsample = 1.
  - The KFF and Health Affairs unwinding measures exist only for four named states (ha_usds_assist); their state tables are not extractable.
  - CMS reporting ends with the July 2024 cohort, so 2024Q4 onward is a true zero in the CMS-based measures.
