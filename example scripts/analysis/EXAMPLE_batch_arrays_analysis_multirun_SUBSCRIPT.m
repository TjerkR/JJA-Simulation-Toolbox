%% Analysis script for batch of arrays - sub script
% Tjerk Reintsema
% 26-11-2021
%
% THIS SCRIPT SHOULD *NOT* BE RAN ON ITS OWN, INSTEAD USE THE MAIN SCRIPT.
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



%% Loading files

if TAKE_AVERAGE || runs == 1
    filestruct = dir(fullfile(datafolder, '*.mat'));
else
    filestruct = dir(fullfile(datafolder, strcat("*_run",string(run),".mat")));
end
files = extractfield(filestruct, 'name');
files = natsort(files);

if length(files) < length(x)
    x = x(1:length(files));
end

%% Calculation f_average between peaks

Ic_array = cell(length(files),1);
zeropeak_locs = zeros(length(files),1);
Ic_max = zeros(length(files),1);
f_average = zeros(length(files),1);
df = zeros(length(files),1);
firstpeak_Ics = zeros(length(files),1);
lastpeak_Ics = zeros(length(files),1);
for n = 1:length(files)
    if TAKE_AVERAGE || runs == 1
        i = ceil(n/runs);
        run = mod(n-1, runs) + 1;
    else
        i = n;
    end
    load(fullfile(datafolder, files{n}))
    if PRINT_NUMBER
        if exist('run', 'var') && runs ~= 0
            fprintf("\nfile: \t%d %d\n\n", x(i), run)
        else
            fprintf("\nfile: \t%d \n\n", x(i))
        end
    end
    
    % Used for errorbars
    df(n) = abs(f_list(end) - f_list(1))/(length(f_list)-1);
    
    % Plot Ic(B) curves and Ic_max(B) for each array
    array_title = erase(files{n}, ".mat");
    if SHOW_ARRAY_PLOTS
        [fig1, ax1] = plot_IcB(Ic_f_list, f_list, nHole_list, 'figno', 1000+n, 'titlestring', strrep(array_title, "_", "\_"));
        [fig2, ax2] = plot_IcB_max(Ic_f_max, f_list, 'figno', 2000+n, 'titlestring', strrep(array_title, "_", "\_"));
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
        Ic_array{n} = visualize_Ic(Ics{i}, 'L', L, 'figno', 3000+n, 'titlestring', strrep(array_title, "_", "\_"));
        
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

        f_average(n) = abs(last_f - first_f) / numel(peaks);
        if PRINT_PEAKS
            fprintf('peaks: \t%.0f\n\n', numel(peaks))
        end

        zeropeak_locs(n) = f_list(locs(ceil(length(nHole_list)/2)));
        
        firstpeak_Ics(n) = Ic_f_max(locs(1));
        lastpeak_Ics(n) = Ic_f_max(locs(end));
    elseif length(locs) == 1
        if f_list(locs(1)) > 0
            firstpeak_Ics(n) = Ic_f_max(locs(1));
            lastpeak_Ics(n) = NaN;
        else
            firstpeak_Ics(n) = NaN;
            lastpeak_Ics(n) = Ic_f_max(locs(1));
        end
        if PRINT_PEAKS
            fprintf('peaks: \t%.0f\n\n', numel(peaks))
        end
    else
        if PRINT_PEAKS
            fprintf('peaks: \t%.0f\n\n', 0)
        end
        firstpeak_Ics(n) = NaN;
        lastpeak_Ics(n) = NaN;
        zeropeak_locs(n) = NaN;
    end
    
    if PRINT_MEAN
        fprintf('mean: \t%.4f\n\n', mean(Ics{i}))
    end
    
    Ic_max(n) = max(Ic_f_max);
    
    if CLOSE_PLOTS
        close all
    end
end

%% EXTRA STEP: AVERAGING
if TAKE_AVERAGE
    f_average = mean(reshape(f_average,runs,[]),1)';
    df = mean(reshape(df,runs,[]),1)' / sqrt(runs);
    Ic_max = mean(reshape(Ic_max,runs,[]),1)';
    zeropeak_locs = mean(reshape(abs(zeropeak_locs),runs,[]),1)'; % taking absolute value now
    firstpeak_Ics = mean(reshape(firstpeak_Ics,runs,[]),1)';
    lastpeak_Ics = mean(reshape(lastpeak_Ics,runs,[]),1)';
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
    clf
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
    clf
    plot(x, zeropeak_locs, '.', 'MarkerSize', 15)
    title(strcat("Zero peak offset as a function of ", xvariable_titlestring))
    xlabel(xlabel_string)
    ylabel("zero peak offset (absolute)")
    if SAVE_OFFSET
        ax = gca;
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "offset_abs_", array_filename_description, ".png"), "Resolution", 400)
    end
end

%% Plot 4: Tilt

if SHOW_FIRSTLAST
    figure(4004)
    clf
    hold on
    plot(x, firstpeak_Ics-mean([firstpeak_Ics, lastpeak_Ics], 2), '.-')
    plot(x, lastpeak_Ics-mean([firstpeak_Ics, lastpeak_Ics], 2), '.-')
    title("Relative I_c of first and last peaks")
    % title(strcat("Relative I_c of first and last peaks as a function of ", xvariable_titlestring))
    xlabel(xlabel_string)
    ylabel("I_c (average of both peaks subtracted)")
    legend(["first peak", "last peak"], 'Location', 'Best')
    if SAVE_FIRSTLAST
        ax = gca;
        if not(isfolder(figuresfolder))
            mkdir(figuresfolder)
        end
        exportgraphics(ax, strcat(figuresfolder, "tilt_rel_", array_filename_description, ".png"), "Resolution", 400)
    end
end
