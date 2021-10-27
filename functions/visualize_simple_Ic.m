function visualize_simple_Ic(array, NameValueArgs)
arguments
    array double
    NameValueArgs.L = 0
    NameValueArgs.figno = 0
    NameValueArgs.titlestring = ""
end

L = NameValueArgs.L;
figno = NameValueArgs.figno;
titlestring = NameValueArgs.titlestring;

N = size(array,2);
rows = N*2;
cols = N;
% removing hole
if L > 0
    hole_edge = (N-L)/2+1;
    array((hole_edge*2-1):(rows+2-hole_edge*2), hole_edge:(cols+1-hole_edge)) = 0;
end

array = flipud(array); % flip or not??
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

end

