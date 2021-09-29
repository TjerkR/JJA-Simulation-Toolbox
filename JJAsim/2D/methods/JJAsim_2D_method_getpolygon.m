function [parea,cx,cy] = JJAsim_2D_method_getpolygon(x,y)
parea = getArea2D(x,y);
[cx,cy] = getCentroid2D(x,y);
end

function area = getArea2D(x,y)
area = sum(x.*circshift(y,[-1,0])-y.*circshift(x,[-1,0]))/2;
end

function [cx,cy] = getCentroid2D(x,y)
area = sum(x.*circshift(y,[-1,0])-y.*circshift(x,[-1,0]))/2;
if area ~= 0
    cx = sum((x + circshift(x,[-1,0])).*(x.*circshift(y,[-1,0])-y.*circshift(x,[-1,0])))/6/area;
    cy = sum((y + circshift(y,[-1,0])).*(x.*circshift(y,[-1,0])-y.*circshift(x,[-1,0])))/6/area;
else
    cx = mean(x,1);
    cy = mean(y,1);
end
end

