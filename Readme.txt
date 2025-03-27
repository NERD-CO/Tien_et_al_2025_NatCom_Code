Description of the publicly available code associated with the manuscript: "Neurons in human motor thalamus encode reach kinematics and positional error related to braking" by Tien et al, submitted to Nature Communications March 2025. Hosted at https://github.com/NERD-CO/Tien_et_al_2025_NatCom_Code

Rex Tien - March 27, 2025
Rex.Tien@cuanschutz.edu

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This code is intended to reproduce all results and figures in the manuscript. The code operates on the publicly available data hosted at: https://osf.io/3unka/?view_only=4e7e0d5cee0e4dd9a567fa637644b3b9

Dependencies: MATLAB with parallel computing toolbox.

This code has been tested with MATLAB R2022a on a system running Windows 11 with 32GB RAM. It should work with any recent MATLAB release on Windows.

Instructions:
To run the code reproducing all results and figures:
1. Download the data available at: https://osf.io/3unka/?view_only=4e7e0d5cee0e4dd9a567fa637644b3b9.
2. Run the main function "Tien_et_al_2025_NatCom_Code.m" with the directory where you downloaded the data as input.
	For example, in MATLAB, run: Tien_et_al_2025_NatCom_Code('C:\Users\username\Downloads\Data')

The computationally intensive processing steps have been pre-executed for speed. By default the main function will load the precalculated data and report the results and produce the plots from the manuscript. In order to re-execute the processing steps, edit the file "Tien_et_al_2025_NatCom_Code.m" line 21 to set "recalculate" to "true."

The processed data fields "StretchBoot," "LagBoot," and "FBLBoot" may not reproduce exactly due to random number generator usage.

The functions contained in the repository are described below.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Main function:

"Tien_et_al_2025_NatCom_Code.m" - main function that loads data, optionally calls functions to regenerate the processed data, reports results and reproduces plots.


Data processing functions:

"get_Stretch.m" - calculates peri-reach firing rates time-stretched to align to both Reach Start and Reach Ends

"get_StretchBoot.m"* - calculates a random distribution of firing rates for significance testing.

"get_StretchSig.m" - calculates significance of modulation of trial-averaged time-stretched peri-reach firing rates.

"get_Realtime.m" - calculates peri-reach and peri-event firing rates with realtime sampling.

"get_RealtimeSig.m" - calculates significance of modulation of trial-averaged realtime peri-reach firing rates.

"get_Directionality.m" - calculates significance of directional tuning of peri-reach firing rates.

"get_RegressionBuffer.m" - generates data buffers for reach regression analyses.

"get_Lag.m" - calculates lagged reach multiple linear regression results.

"get_LagBoot.m"* - calculates a random distribution of reach multiple linear regression results for significance testing.

"get_LagSig.m" - calculates significance of reach regression results.

"get_ShapLag.m"* - performs Shapley Value Decomposition for reach regressions.

"get_FBRegressionBuffer.m" - generates data buffers for start and end windowed regression analyses.

"get_FBL.m" - calculates lagged start and end windowed multiple linear regression results.

"get_FBLBoot.m"* - calculates a random distribution of start and end windowed multiple linear regression results for significance testing.

"get_FBLSig.m" - calculates significance of start and end windowed regression results.

* These functions are compute-intensive and take a long time to run on most systems.


Plotting functions:

"plot_speed_segment.m" - plots a fingertip speed segment (as in Figure 1a).

"plot_3DKin.m" - plots the 3D fingertip position from a session (as in Figure 1e).

"plot_recording_locations.m" - plots the recording locations (track and depth) of all analyzed units (Figure 2a).

"plot_Raster.m" - plots a unit's peri-reach spike raster colored by reach direction (as in Figure 2b,c).

"plot_DirMeanFRs.m" - plots a unit's peri-reach firing rates aligned to peak speed, colored by reach direction (as in Figure 2d,e).

"plot_report_Stretch.m" - plots peri-reach trial-averaged FRs and peri-reach modulation significance results (Figure 3). Also plots recording locations of significant peri-reach modulated units (Extended Data Figure 5a,b).

"plot_PCRF.m" - plots a unit's spatial activation map (as in Figure 4a).

"plot_report_Directionality.m" - plots significant peri-reach directional tuning results (Figure 4b,c). Also plots recording locations of significant peri-reach directionally tuned units (Extended Data Figure 5c,d).

"plot_report_ShapLag.m" - plots regression results (Figure 5). Also plots recording locations of significant kinematic encoding units (Extended Data Figure 5e,f).

"plot_Realtime.m" - plots peri-reach and peri-event trial-averaged realtime FRs and realtime peri-reach modulation significance results (Extended Data Figures 1 and 2).

"plot_FBL.m" - plots start and end windowed regression results (Extended Data Figure 3).

"plot_Cutoffs_Reach.m" - plots the prevalence of peri-reach modulation significance for different required significance durations (Extended Data Figure 4a).

"plot_Cutoffs_Directionality.m" - plots the prevalence of peri-reach directional tuning significance for different required significance durations (Extended Data Figure 4b).


Utility functions:

"GaussSmooth_Arbitrary.m" - smooths and resamples by convolving with a Gaussian kernel at desired sample points.

"make_save_struct.m" - places variables into a structure for saving.

"my_normalize.m" - does z-score normalization while accounting for NaN values.

"Spikes2FIR_Arbitrary.m" - calculates fractional interval firing rates at desired sample points.

"target_idx_swap.m" - re-indexes targets for more intuitive plotting.



