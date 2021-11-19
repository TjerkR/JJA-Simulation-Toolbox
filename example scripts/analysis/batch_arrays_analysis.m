%% Plot average distance between peaks for batch of arrays
% Tjerk Reintsema
% 19-11-2021
%
% Analysis script for Ic(B) data for a batch of different arrays.
% Functionality includes plotting individual and maximum Ic(B) curves for
% each array, showing the distribution of the critical currents for each
% array and plotting the average distance between peaks, maximum critical
% current, offset and first and last peak critical currents as a function 
% of an array parameter.
% -> adapted from:
% "\[2021-11-05] 28a20h 0fracN...\" on 19-11-2021
close all
clearvars

%% INPUT / OPTIONS

%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel_string = "I_c imbalance divide position";
xvariable_titlestring = "array I_c imbalance divide position";
array_filename_description = "28a20h_fracN_al_lohi_varin0p001"; 

fracNs = 0:1:28;
x = fracNs; 
L = 20;

figuresfolder = '.\figures\';
visualizationsfolder = '.\visualizations\';
datafolder = '.\data\';

%%%% PLOTTING / SAVING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SHOW_ALL = false;
SAVE_ALL = false;

% Array plots
SHOW_ARRAY_PLOTS = SHOW_ALL;
SAVE_ARRAY_PLOTS = SAVE_ALL;

% Array visualizations
SHOW_VISUALIZATIONS = SHOW_ALL;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading files
filestruct = dir(fullfile(datafolder, '*.mat'));
files = extractfield(filestruct, 'name');
files = natsort(files);
x = x(1:length(files));

%% Calculation f_average between peaks

Ic_array = cell(length(files),1);
zeropeak_locs = zeros(length(files),1);
Ic_max = zeros(length(files),1);
f_average = zeros(length(files),1);
df = zeros(length(files),1);
firstpeak_Ics = zeros(length(files),1);
lastpeak_Ics = zeros(length(files),1);
for i = 1:length(files)
    load(fullfile(datafolder, files{i}))
    disp(x(i))
    
    % Used for errorbars
    df(i) = abs(f_list(end) - f_list(1))/(length(f_list)-1);
    
    % Plot Ic(B) curves and Ic_max(B) for each array
    array_title = erase(files{i}, ".mat");
    if SHOW_ARRAY_PLOTS
        [fig1, ax1] = plot_IcB(Ic_f_list, f_list, nHole_list, 'figno', 1000+i, 'titlestring', strrep(array_title, "_", "\_"));
        [fig2, ax2] = plot_IcB_max(Ic_f_max, f_list, 'figno', 2000+i, 'titlestring', strrep(array_title, "_", "\_"));
        IcB_xlim = ax2.XLim;
        IcB_ylim = ax2.YLim;
        ax1.XLim = IcB_xlim;
        ax1.YLim = IcB_ylim;

        if SAVE_ARRAY_PLOTS
            if not(isfolder(figuresfolder))
                mkdir(figuresfolder)
            end
            exportgraphics(ax1, strcat(figuresfolder, array_title, "_all.png"), "Resolution", 400)
            exportgraphics(ax2, strcat(figuresfolder, array_title, "_max.png"), "Resolution", 400)
        end
    end
        
    % Show visualization of critical currents
    if SHOW_VISUALIZATIONS
        Ic_array{i} = visualize_Ic(Ics{i}, 'L', L, 'figno', 3000+i, 'titlestring', strrep(array_title, "_", "\_"));
        fprintf('mean: %.4f\n\n', mean(Ics{i}))
        
        if SAVE_VISUALIZATIONS
            if not(isfolder(visualizationsfolder))
                mkdir(visualizationsfolder)
            end
            ax = gca;
            exportgraphics(ax, strcat(visualizationsfolder, array_title, "_visual.png"), "Resolution", 400)
        end
    end
        
    % Find peaks and calculate f_average
    [peaks, locs] = findpeaks(Ic_f_max, 'MinPeakHeight', max(Ic_f_max)/4);
    if SHOW_ARRAY_PLOTS
        for loc = locs
            xline(f_list(loc))
        end
    end
    
    if ~isempty(locs)
        first_f = f_list(locs(1));
        last_f = f_list(locs(end));

        f_average(i) = abs(last_f - first_f) / numel(peaks);
        fprintf('peaks: %.0f\n\n\n', numel(peaks))

        zeropeak_locs(i) = f_list(locs(ceil(length(nHole_list)/2)));
        
        firstpeak_Ics(i) = Ic_f_max(locs(1));
        lastpeak_Ics(i) = Ic_f_max(locs(end));
    elseif length(locs) == 1
        if f_list(locs(1)) > 0
            firstpeak_Ics(i) = Ic_f_max(locs(1));
            lastpeak_Ics(i) = NaN;
        else
            firstpeak_Ics(i) = NaN;
            lastpeak_Ics(i) = Ic_f_max(locs(1));
        end
    else
        fprintf('peaks: %.0f\n\n\n', 0)
        firstpeak_Ics(i) = NaN;
        lastpeak_Ics(i) = NaN;
        zeropeak_locs(i) = NaN;
    end
    
    Ic_max(i) = max(Ic_f_max);
end


%% Plot 1: f_average

if SHOW_F_AVERAGE
    titlestring = strcat("f_{average} between peaks as a function of ", xvariable_titlestring);
    xlims = 'auto';
    ylims = 'auto';
    figno = 4001;

    figure(figno)
    clf

    if ERRORBARS
        errorbar(x, f_average, df, '.', 'MarkerSize', 15)
    else
        plot(x, f_average, '.', 'MarkerSize', 15)
    end
    title(titlestring)
    xlabel(xlabel_string)
    ylabel("f_{average}")

    ax = gca;
    ax.YRuler.Exponent = 0;
    ax.XTickMode = 'manual';
    ax.XLimMode = 'manual';
    box on

    xlim(xlims)
    ylim(ylims)

    if SAVE_F_AVERAGE
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "f_average_", array_filename_description, ".png"), "Resolution", 400)
    end
end

%% Plot 2: Ic_max

if SHOW_IC_MAX
    figure(4002)
    plot(x, Ic_max, '.', 'MarkerSize', 15)
    title(strcat("Maximum critical current as a function of ", xvariable_titlestring))
    xlabel(xlabel_string)
    ylabel("maximum I_c")
    if SAVE_IC_MAX
        ax = gca;
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "Ic_max_", array_filename_description, ".png"), "Resolution", 400)
    end
end

%% Plot 3: Offset

if SHOW_OFFSET
    figure(4003)
    plot(x, zeropeak_locs, '.', 'MarkerSize', 15)
    title(strcat("Zero peak offset as a function of ", xvariable_titlestring))
    xlabel(xlabel_string)
    ylabel("zero peak offset")
    if SAVE_OFFSET
        ax = gca;
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "offset_", array_filename_description, ".png"), "Resolution", 400)
    end
end

%% Plot 4: First-last peak

if SHOW_FIRSTLAST
    figure(4004)
    hold on
    plot(x, firstpeak_Ics, '.-')
    plot(x, lastpeak_Ics, '.-')
    title(strcat("I_c of first and last peaks as a function of ", xvariable_titlestring))
    xlabel(xlabel_string)
    ylabel("I_c of first and last peak")
    legend(["first peak", "last peak"], 'Location', 'Best')
    if SAVE_FIRSTLAST
        ax = gca;
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "tilt_", array_filename_description, ".png"), "Resolution", 400)
    end
end
