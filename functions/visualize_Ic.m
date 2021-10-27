function array = visualize_Ic(Ic, NameValueArgs)
arguments
    Ic (:,1) double
    NameValueArgs.L = 0
    NameValueArgs.figno = 0
    NameValueArgs.titlestring = ""
end
L = NameValueArgs.L;
figno = NameValueArgs.figno;
titlestring = NameValueArgs.titlestring;


array = flipud(Ic_list_to_array(Ic, L)); % flip or not??
array = [array array(:,end)];
array = [array; array(end,:)];

if NameValueArgs.figno == 0
    figure;
else
    figure(figno);
end
clf

pcolor(array)
title(titlestring)

array = Ic_list_to_array(Ic, L);
end