%% Example - Plot average distance between peaks for batch of arrays
% Tjerk Reintsema
% 18-10-2021
%
% Plots and optionally saves the average distance between peaks as a
% function of Ic variance.
% (Taken from batch_arrays_analysis_149285)
close all
clearvars

SAVE_ARRAY_PLOTS = true;
SAVE_OVERVIEW_PLOT = true;
ERRORBARS = true;

folder = '.\data\';
variances = linspace(0, 0.10, 21);

filestruct = dir(fullfile(folder, '*.mat'));
files = extractfield(filestruct, 'name');
files = natsort(files);

variances = variances(1:length(files));

%% Calculation f_average between peaks

Ic_max = zeros(length(files),1);
f_average = zeros(length(files),1);
df = zeros(length(files),1);
for i = 1:length(files)
    load(fullfile(folder, files{i}))
    
    % Used for errorbars
    df(i) = abs(f_list(end) - f_list(1))/(length(f_list)-1);
    
    % Plot Ic_max(B) for each array
    array_title = erase(files{i}, ".mat");
    [fig1, ax1] = plot_IcB(Ic_f_list, f_list, nHole_list, 'figno', i, 'titlestring', strrep(array_title, "_", "\_"));
    [fig2, ax2] = plot_IcB_max(Ic_f_max, f_list, 'figno', 100+i, 'titlestring', strrep(array_title, "_", "\_"));
    
    if SAVE_ARRAY_PLOTS
        if not(isfolder(".\figures\"))
            mkdir(".\figures\")
        end
        exportgraphics(ax1, strcat(".\figures\", array_title, "_all.png"), "Resolution", 400)
        exportgraphics(ax2, strcat(".\figures\", array_title, "_max.png"), "Resolution", 400)
    end
    
    % Print number of peaks for each array
%     disp(num_peaks(Ic_f_max))
    
    % Find peaks and calculate f_average
    % OPTION 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Just find all peaks %%%%%%%%%%%%%%%%%%%%%%%%%%
%     [peaks, locs] = findpeaks(Ic_f_max, 'MinPeakHeight', max(Ic_f_max)/4);
%     for loc = locs
%         xline(f_list(loc))
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % OPTION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Try to get only "normal" peaks%%%%%%%%%%%%%%%%
    location = ceil(16/30* length(Ic_f_max));
    [peaks, locs] = findpeaks(Ic_f_max(location:end), 'MinPeakHeight', max(Ic_f_max)/4);
    for loc = locs
        xline(f_list(location+loc))
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    first_f = f_list(locs(1));
    last_f = f_list(locs(end));

    f_average(i) = abs(last_f - first_f) / numel(peaks);
    disp(numel(peaks))
    Ic_max(i) = max(Ic_f_max);
end



%% Plotting

filestring = "f_average 34a30h square disordered";
titlestring = "f_{average} between peaks for different Ic disorder variances";
xlims = 'auto';
ylims = 'auto';
figno = 1111; % use number >1000 to prevent overwriting

figure(figno)
clf

if ERRORBARS
    errorbar(variances, f_average, df, '.', 'MarkerSize', 15)
else
    plot(variances, f_average, '.', 'MarkerSize', 15)
end
title(titlestring)
xlabel("I_c variance")
ylabel("f_{average}")

ax = gca;
ax.YRuler.Exponent = 0;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

xlim(xlims)
ylim(ylims)

if SAVE_OVERVIEW_PLOT
    if not(isfolder(".\figures\"))
        mkdir(".\figures\")
    end
    exportgraphics(ax, strcat(".\figures\", filestring, ".png"), "Resolution", 400)
end