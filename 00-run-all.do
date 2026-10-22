*==============================================================================
* 00-run-all.do
*
* Runs the whole package. Edit the paths and switches here once; the numbered files read these globals.
*
* Order
*   (a) controls\build_politics_panel.do, build_unwinding_panel.do, merge_state_controls.do
*       (the state controls package: NCSL politics, CMS unwinding, Health Affairs flag) -> controls_state_yq.dta
*   01  01-build-state-controls.do   merges everything onto the state-quarter file and prints the VIFs
*   02  02-agg-util-did.do           aggregate utilization DiD (4 outcomes x S1-S4 x main/ext)
*   03  03-ind-did.do                individual imputation DiD: bmi | medicaid == 1 and medicaid100 (S1-S4 x main/ext)
*   04  04-ind-iv.do                 individual IV: glp instruments obesity and placebo utilization and spending
*
* Each of 02-04 writes its own log with every specification and the comparison tables:
*   agg_util_did24.log  agg_util_did25.log  ind_did24.log  ind_did25.log  ind_iv24.log  ind_iv25.log
* (24 = main sample, 25 = extended sample through 2025)
*
* Packages: ssc install did_imputation, replace   (did_imputation needs reghdfe and ftools: ssc install reghdfe ; ssc install ftools)
* Not run: Stata was not available when this was written.
*==============================================================================
clear all
set more off

*---- paths --------------------------------------------------------------------
global root "C:/Users/jsemprini/OneDrive - Des Moines University/9-Internal-Collab/a-Devine_GLP1"
global code "C:\Users\jsemprini\OneDrive - Des Moines University\9-Internal-Collab\a-Devine_GLP1\7-package\estimation"            // folder that holds these do files

*---- switches -------------------------------------------------------------------
global runs      main ext                  // main = 2021Q1-2024Q4 ; ext = through 2025Q4
global specs     S1 S2 S3 S4               // S1 unadjusted, S2 + unwinding, S3 + enrollment, S4 + policy controls ALL
global fixsample 1                         // 1 = same observations in every specification (recommended for comparing)
global outsel    bmi medicaid100           // individual outcomes for 03 and 04
global dyn       1                         // 03 only: 1 = also run the horizons / pretrends models (slow on the medicaid sample)

do "$code/01-build-state-controls.do"
do "$code/02-agg-util-did.do"
do "$code/03-ind-did.do"
do "$code/04-ind-iv.do"
