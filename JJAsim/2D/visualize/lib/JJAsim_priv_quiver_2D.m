function [p,xap,yap] = JJAsim_priv_quiver_2D(ah,U,V,X,Y,scale,stretch,wscale,plotQ)
if nargin < 6
    scale = 1;
end
if nargin < 7
    stretch = 1;
end
if nargin < 8
    wscale = 1;
end
if nargin < 9
    plotQ = true;
end
x = (X + U/2)';
y = (Y + V/2)';
L = sqrt(U.^2+V.^2)'*stretch;
alpha = atan2(V,U)';
W = ((L/max(L)).^0.5)*max(L)*wscale;
pl = 0.7-0.3*L/max(L);
[xap,yap]= parrow(x,y,W*0.2*scale,W*0.5*scale,L.*pl*scale,L,alpha);
if plotQ
    p = patch(xap,yap,[0,0,1]);
    p.EdgeColor = [0,0,1];
    p.Parent = ah;
else
    p = [];
end
end


function [xap,yap]= parrow(x,y,w,W,l,L,alpha)
xa = [-L/2;L/2-l;L/2-l;L/2;L/2-l;L/2-l;-L/2];
ya = [w;w;W;zeros(size(W));-W;-w;-w]/2;
xap = bsxfun(@times,cos(alpha),xa)-bsxfun(@times,sin(alpha),ya);
yap = bsxfun(@times,cos(alpha),ya)+bsxfun(@times,sin(alpha),xa);
xap = bsxfun(@plus,xap,x);
yap = bsxfun(@plus,yap,y);
end