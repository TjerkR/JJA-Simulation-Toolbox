%% Analysis script for batch of arrays - main script
% Tjerk Reintsema
% 26-11-2021
%
% THIS IS THE MAIN SCRIPT THAT NEEDS TE BE RAN.
%
% Analysis script for Ic(B) data for a batch of different arrays.
% Functionality includes plotting individual and maximum Ic(B) curves for
% each array, showing the distribution of the critical currents for each
% array and plotting the average distance between peaks, maximum critical
% current, offset and first and last peak critical currents as a function 
% of an array parameter.
% Includes option for averaging multiple runs.
%
% -> adapted from:
% "\[2021-11-18] 28a20h ... multiple runs\" on 22-11-2021
% close all
clearvars

%% INPUT / OPTIONS

%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel_string = "I_c imbalance divide position";
xvariable_titlestring = "array I_c imbalance divide position";
array_filename_description = "28a20h_fracN_al_lohi_nodisorder"; 

fracNs = 0:1:28;
x = fracNs; 
L = 20;

figuresfolder = '.\figures EWI\';
visualizationsfolder = '.\visualizations EWI\';
datafolder = '.\data EWI\';

runs = 1; % always needs to be set!
TAKE_AVERAGE = false;
run = 0;

%%%% PLOTTING / SAVING / PRINTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHOW_ALL = true;
SAVE_ALL = false;

% Array plots
SHOW_ARRAY_PLOTS = false;
SAVE_ARRAY_PLOTS = SAVE_ALL;
CLOSE_PLOTS = false; % enable this to close array plots after saving

% Array visualizations
SHOW_VISUALIZATIONS = false;
SAVE_VISUALIZATIONS = SAVE_ALL;

% Overview plots
SHOW_F_AVERAGE = SHOW_ALL;
SAVE_F_AVERAGE = SAVE_ALL;
ERRORBARS = true;

SHOW_IC_MAX = SHOW_ALL;
SAVE_IC_MAX = SAVE_ALL;

SHOW_OFFSET = SHOW_ALL;
SAVE_OFFSET = SAVE_ALL;

SHOW_FIRSTLAST = SHOW_ALL;
SAVE_FIRSTLAST = SAVE_ALL;

PRINT_NUMBER = true;
PRINT_MEAN = true;
PRINT_PEAKS = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculations and plotting

batch_arrays_analysis_multirun_SUBSCRIPT;
