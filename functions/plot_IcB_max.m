function [fig, ax] = plot_IcB_max(Ic_f_max, f_list, NameValueArgs)
%UNTITLED Summary of this function goes here
%
arguments
    Ic_f_max (1,:) double
    f_list (1,:) double
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

plot(f_list, Ic_f_max, '.')
xlabel('f')
ylabel('I_c')

title(titlestring)
xlim(xlims)
ylim(ylims)

ax = gca;
ax.XTickMode = 'manual';
ax.XLimMode = 'manual';
box on

end

