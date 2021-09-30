function [fig, ax] = plot_IcB_max(Ic_f_max, f_list, figno, titlestring, xlims, ylims)
%UNTITLED Summary of this function goes here
%

figure(figno)
clf

plot(f_list, Ic_f_max, '.')

fig = figure(figno);
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

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

% exportgraphics(ax, 'Title.png', 'Resolution', 400)

end

