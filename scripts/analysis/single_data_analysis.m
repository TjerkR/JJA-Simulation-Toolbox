close all
clearvars

% file = 'C:\Users\tjerk\Documents\MATLAB\Graduation Project\Reproduction thesis Julia\data analysis\FOLDER\file.mat';
file = '.\data analysis\single geometries\20a6h_v2.mat';

load(file);

%%

for i = 1:size(Ic_f_list, 1)
    figure(10);
    hold on
    plot(f_list, Ic_f_list(i, :),'.')
end

fig1 = figure(10);
legend(strcat({'v = '},string(num2cell(nHole_list))))
title('I_c of a 20 by 20 array, with a 6 by 6 hole')
xlabel('f')
ylabel('I_c')
xlim([-0.08, 0.08])
ylim([0, 0.7])

legend('boxoff')
legend('FontSize', 6)

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

xl = xlim;
yl = ylim;
xt = xticks;
yt = yticks;

exportgraphics(ax, '20a6h_all.png', 'Resolution', 400)



Ic_f_max = max(Ic_f_list);
fig2 = figure(20);
plot(f_list, Ic_f_max, '.')
xlim(xl)
ylim(yl)
xticks(xt)
yticks(yt)
title('I_c of a 20 by 20 array, with a 6 by 6 hole')
xlabel('f')
ylabel('I_c')

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

exportgraphics(ax, '20a6h_max.png', 'Resolution', 400)
