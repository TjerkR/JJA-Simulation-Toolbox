%% Example - Plot average distance between peaks for batch of arrays
% Tjerk Reintsema
% 01-10-2021
%
% Plots and optionally saves the average distance between peaks as a
% function of hole size.
close all
clearvars

folder = '.\data\';
holes = 8:1:13;

filestruct = dir(fullfile(folder, '*.mat'));
files = extractfield(filestruct, 'name');
files = natsort(files);

%% Calculation f_average between peaks

f_average = zeros(length(files),1);
df = zeros(length(files),1);
for i = 1:length(files)
    load(fullfile(folder, files{i}))
    
    % Used for errorbars
    df(i) = abs(f_list(end) - f_list(1))/(length(f_list)-1);
    
    % Plot Ic_max(B) for each array
    plot_IcB_max(Ic_f_max, f_list, 'figno', i, 'titlestring', strrep(files{1}, "_", "\_"))
    
    % Print number of peaks for each array
    disp(num_peaks(Ic_f_max))
    
    % Find peaks and calculate f_average
    [peaks, locs] = findpeaks(Ic_f_max, 'MinPeakHeight', max(Ic_f_max)/4);
    first_f = f_list(locs(1));
    last_f = f_list(locs(end));

    f_average(i) = abs(last_f - first_f) / numel(peaks);
end

%% Plotting

SAVE_FIGURES = true;
ERRORBARS = true;
filestring = "f_average (L+2)a8-13h square";
titlestring = "f_{average} between peaks for different hole sizes, N = L+2";
xlims = 'auto';
ylims = 'auto';
figno = 111; % use high number to prevent overwriting

figure(figno)
clf

if ERRORBARS
    errorbar(holes, f_average, df, '.', 'MarkerSize', 15)
else
    plot(holes, f_average, '.', 'MarkerSize', 15)
end
title(titlestring)
xlabel("hole size")
ylabel("f_{average}")

ax = gca;
ax.YRuler.Exponent = 0;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

xlim(xlims)
ylim(ylims)