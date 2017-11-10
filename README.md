# MODGARNER
This project contains data and code for the modified Garner experiments published in:
Lin & Little (2017)

Additional code is included for anlaysis of published in:
Little, D. R., Wang, T. & Nosofsky, R. (2016). Sequence-sensitive exemplar and decision-bound accounts of speeded-classification performance in a modified Garner-tasks paradigm. Cognitive Psychology, 89, 1-38. doi: 10.1016/j.cogpsych.2016.07.001

Archive created 10-Nov-17

========================================================================
Folders:
========================================================================
Experimental Code
------------------------------------------------------------------------
ModGarnerSeparable
- Contains MATLAB code for running modified Garner task with Brightness and Saturation in different patches
-- RUN_SEPARABLE.m is the main script
-- Data was initial collected as part of honours project for Deborah Lin

MODIFIED GARNER BOXCAR/Experiment
- Contains MATLAB code for running modified Garner task with Line position and Saturation
-- RUN_SEPARABLE_BOXCAR.m is the main script

========================================================================
Data
- Contains Raw Data and analysis code from Lin2017
------------------------------------------------------------------------
ModGarnerSeparable
- \rawdata\ contains the raw data from the ModGarnerSeparable Experiment
-- columns are: 
cols = {'sub', 'con', 'side', 'sess', 'task', 'prac', 'trl', 'itm', 'leftBri', 'leftSat', 'rightBri', 'rightSat', 'corrResp', 'resp', 'corr', 'rt', 'exc'};
-- side indicates the side of the relevant dimension (left or right)
-- task [1 control A; 2 control B; 3 correlated A; 4 correlated B; 5 filtering; 6 filtering]

- \analysis\ contains the code to plot the sequential effects data and to pre-process the data for analysis of the sequential effects

- \data\ contains the data after pre-processing
-- columns are:
cols = {'sub', 'task', 'sess', 'trl', 'item', 'acc', 'rts', 'exc', 'igp'}
-- exc = 1 if the trial should be excluded
-- igp = 1 if the trial is the first trial in a block (and, hence, no preceeding trial)

- Subject numbers correspond to the paper as follows:
101-104 = B1-B4, 201-204 = S1-S4

- \Sequential Effects Analysis - Multilevel Regression Analysis\ contains the data and [R] code to perform the Bayesian Multilevel Regression analysis of the sequential effects
- \Sequential Effects Analysis - Plot\ contains the code to plot the sequential effects
- \Individual RT analysis\ contains the JAGS code to apply the log RT model to the individual RTs

ModGarnerBoxcar
- \rawdata\ as above
-- columns are:
cols = {'sub', 'con', 'sess', 'task', 'prac', 'trl', 'itm', 'linepos', 'saturation',  'corrResp', 'resp', 'corr', 'rt', 'exc'};

- \analysis\ as above
- \data\ as above
- Subject numbers correspond to the paper as follows:
101-104 = L1-L4, 205, 202-204 = S1-S4
- \Sequential Effects Analysis - Multilevel Regression Analysis\ as above
- \Sequential Effects Analysis - Plot\ as above
- \Individual RT analysis\ as above

ModGarnerIntegral
- \rawdata\ as above
-- columns are:
cols = {'sub', 'con', 'sess', 'task', 'prac', 'trl', 'itm', 'saturation', 'brightness',  'corrResp', 'resp', 'corr', 'rt', 'exc'};

- \data\ as above

AverageOverallRTanalysis
- contains JAGS code to analyse the overall RTs for the Integral, Separable & Boxcar Experiments


========================================================================
Modeling
------------------------------------------------------------------------
GARNER SEQ EB-LBA DEMCMC
- Contains code for fitting sequence-sensitive exemplar model using DEMCMC
- Main file is fitDEMCMC_GarnerWithin_v3b.m
-- This will run the DEMCMC estimation
- Use plotPosteriors to view chains
- Use plotPosteriorPredictives7.m to make plots in the supplementary material (need to put fit file in fit folder first)
- Use plotPosteriorsViolinExp to make parameter plots in the supplement
- When fitting integral dimensions, need to change R_METRIC = 2 and SHEPARD_EXPONENT =2 in seqEBLBA.m. Set both = 1 for separable dimensions.
