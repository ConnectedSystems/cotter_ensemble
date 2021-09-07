IHACRES State-based Ensemble - Cotter

[![DOI](https://zenodo.org/badge/388680035.svg)](https://zenodo.org/badge/latestdoi/388680035)


Code and analyses for investigations into a state-based IHACRES ensemble modeling approach featuring time-varying parameter sets.

Written for Julia v1.5.3 using [Streamfall.jl](https://github.com/ConnectedSystems/Streamfall.jl) (v0.19) and [ihacres_nim](https://github.com/ConnectedSystems/ihacres_nim) (v0.4)

All instructions below is run from the `src` directory.

To start:

```bash
# Install necessary packages
$ julia --project=..
julia> instantiate
```

Run files in order:

```bash
$ julia --project=..
julia> include("1_analyze_data.jl")
julia> include("2_sensitivity_analysis.jl")
# ... etc ...
```

---



## Overview of files in `src` directory

| File                            	| Description                                                             	|
|---------------------------------	|-------------------------------------------------------------------------	|
| _common.jl                      	| Defines common configurations, variables, helper functions and packages 	|
| _calib_funcs.jl                 	| Defines calibration functions                                           	|
| 1a_analyze_data.jl              	| Analyze climate data and produce figures                                	|
| 2_sensitivity_analysis.jl       	| Preliminary sensitivity analysis                                        	|
| 3a_calib_baseline_IHACRES.jl    	| Generate baseline instances                                             	|
| 4_analyze_baselines.jl          	| Analysis and plots of baseline instances                                	|
| 5a_cmd_state_10_year_burn_in.jl 	| Generate state-based instances                                          	|
| 6a_build_ensembles_.jl          	| Create ensemble pairs                                                   	|
| 6b_compare_ensembles.jl         	| Comparison of ensemble pairs and plots for paper                        	|


## Overview of files in `data` directory

| File                                     	| Description                                                                                          	|
|------------------------------------------	|------------------------------------------------------------------------------------------------------	|
| 410730_IHACRES.yml                       	| Streamfall network specification file                                                                	|
| sa_results_w_dummy.xlsx                  	| Excel file containing preliminary sensitivity analyses                                               	|
| CAMELS-AUS_410730.csv                    	| Data extracted from CAMELS-AUS dataset                                                               	|
| baselines/*                              	| Streamfall network specification files for calibrated baseline instances                             	|
| calib_params/*                           	| Raw dump of parameter values for ensemble instances                                                  	|
| ensemble_results/*_mix.csv               	| Dump of weight values used for wet, usual, dry states                                                	|
| ensemble_results/*_outflows.csv          	| Outflows estimated by each instance                                                                  	|
| ensemble_results/*_NmKGE_results.csv     	| Overview of performance metrics for each instance                                                    	|
| ensemble_results/state_active_params.csv 	| Time series indicating which state was active for a given time step (for state-based instances only) 	|
| ensemble_results/state_state_mix.csv     	| Weight values for State-State instances                                                              	|


## Overview of `figure` directory

| File                   	| Description                                          	|
|------------------------	|------------------------------------------------------	|
| calib_baseline/*       	| Indicative plots of calibrated baseline instances    	|
| calib_state/*          	| Indicative plots of calibrated state-based instances 	|
| climate/*              	| Plots for climate data analyses                      	|
| ensemble_results/cmd/* 	| Figures used for paper                               	|