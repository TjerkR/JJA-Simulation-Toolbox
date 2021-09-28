close all
clearvars

folder = 'C:\Users\tjerk\Documents\MATLAB\Graduation Project\server results\slurm-130984';
sizes = 52:2:66;

filestruct = dir(fullfile(folder, '*.mat'));
files = extractfield(filestruct, 'name');
files = natsort(files);
     


%%

f_average = zeros(length(files),1);
df = zeros(length(files),1);
for i = 1:length(files)
    load(fullfile(folder, files{i}))
    
    df(i) = abs(f_list(end) - f_list(1))/(length(f_list)-1);
    
    figure(i)
    plot(f_list, Ic_f_max)
    xlabel('f')
    ylabel('I_c')
    title(files{i})
    
    disp(num_peaks(Ic_f_max))
    
    [peaks, locs] = findpeaks(Ic_f_max, 'MinPeakHeight', max(Ic_f_max)/4);
    first_f = f_list(locs(1));
    last_f = f_list(locs(end));

    f_average(i) = abs(last_f - first_f) / numel(peaks);
end

%%
figure(111)
clf
errorbar(sizes, f_average, df, '.', 'MarkerSize', 15)
% plot(sizes, f_average, '.', 'MarkerSize', 15)

ax = gca;
ax.YRuler.Exponent = 0;
% ytickformat('%.0g')
title('f_{average} between peaks for different array sizes (constant 50x50 hole)')
% title('50a array')

% load('C:\Users\tjerk\Documents\MATLAB\Graduation Project\Reproduction thesis Julia\data analysis\f_average(holesize)\fitting\fittedmodel2.mat')
% hold on
% plot(fittedmodel)

xlabel('hole size')
ylabel('f_{average}')

% ylim([0.00039, 0.00054])

% legend('data', 'fitted curve')
% legend('boxoff')
% legend('FontSize', 6)

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

exportgraphics(ax, 'f_average 52-66a50h.png', 'Resolution', 400)
