function [fig, ax] = plot_IcB(Ic_f_list, f_list, nHole_list, figno, titlestring, xlims, ylims)
%UNTITLED Summary of this function goes here
%   TODO: Also plot IcB_max

figure(figno)
clf

for i = 1:size(Ic_f_list, 1)
%     figure(figno);
    hold on
    plot(f_list, Ic_f_list(i, :), '.')
end

fig = figure(figno);
legend(strcat({'v = '}, string(num2cell(nHole_list))))
xlabel('f')
ylabel('I_c')

try
    if ~strcmp(titlestring, "none")
        title(titlestring)
    end
end

try
    if class(xlims) == "double"
        xlim(xlims)
    end
end

try
    if class(ylims) == "double"
        ylim(ylims)
    end
end

legend('boxoff')
legend('FontSize', 6)

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

% exportgraphics(ax, 'Title.png', 'Resolution', 400)

end

