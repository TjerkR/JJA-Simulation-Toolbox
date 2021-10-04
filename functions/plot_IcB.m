function [fig, ax] = plot_IcB(Ic_f_list, f_list, nHole_list, NameValueArgs)
%UNTITLED Summary of this function goes here
%
arguments
    Ic_f_list double
    f_list (1,:) double
    nHole_list (1,:) double
    NameValueArgs.figno = ""
    NameValueArgs.titlestring = ""
    NameValueArgs.xlims = 'auto'
    NameValueArgs.ylims = 'auto'
end

titlestring = NameValueArgs.titlestring;
xlims = NameValueArgs.xlims;
ylims = NameValueArgs.ylims;

if class(NameValueArgs.figno) == "string"
    fig = figure;
else
    fig = figure(NameValueArgs.figno);
end
clf

for i = 1:size(Ic_f_list, 1)
    hold on
    plot(f_list, Ic_f_list(i, :), '.')
end

legend(strcat({'v = '}, string(num2cell(nHole_list))))
xlabel('f')
ylabel('I_c')

title(titlestring)
xlim(xlims)
ylim(ylims)

legend('boxoff')
legend('FontSize', 6)

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

end

