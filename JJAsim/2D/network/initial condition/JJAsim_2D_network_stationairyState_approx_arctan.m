function [th,I,n,nInput] = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f)
%[th,I,n] = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f)
%
%DESCRIPTION
% - Generate an approximate stationairy state for a 2D network with the 
%   arctan approximation, which can be used as an initial guess for 
%   JJAsim_2D_network_stationairyState or initial condition for 
%   JJAsim_2D_network_simulate. 
% - The resulting th will encode the specified vortex configuration.
% - A vortex configuration is specified with a list of coordinate and 
%   vorticity pairs, (x0,y0,n0).
% - The output is in the phase zone z = 0. One can use the function
%   JJAsim_2D_network_method_changePhaseZone to convert it to other 
%   phase zones. 
% - Also returns nInput, which is the vector representation of the input
%   vortex configuration. This need not be equal to n, the vortex
%   configuration of the output th. If they differ, a warning is given.
%
%FIXED INPUT
% array      struct      information about Josephson junction array
% x0         L by W*     x coordinates of the L sites that contain vortices
% y0         L by W*     y coordinates of the L sites that contain vortices
% n0         L by W*     vorticity of the L sites that contain vortices.
%                        Must be integers.
% f          Np* by W*   frustration factor 
%
%OUTPUT
% th         Nj by W     approximate gauge invariant phase difference
% I          Nj by W     approximate junction current
% n          Np by W     vortex configuration of approximate initial condition
% nInput     Np by W     vector representation of input vortex configuration

if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

if mean(mean(abs(round(n0)-n0))) > 1E-10
    error('n0 must be integers');
end

%get sizes
Np = array.Np;
L = max([size(x0,1),size(y0,1),size(n0,1)]);
W = max([size(x0,2),size(y0,2),size(n0,2),size(f,2)]);

%check input
x0 = JJAsim_method_checkInput(x0,'double',[L,W],[0,1],'x0');
y0 = JJAsim_method_checkInput(y0,'double',[L,W],[0,1],'y0');
n0 = JJAsim_method_checkInput(n0,'double',[L,W],[0,1],'n0');
f = JJAsim_method_checkInput(f,'double',[Np,W],[1,1],'f');
if size(x0,2) == 1
    x0 = repmat(x0,1,W);
end
if size(y0,2) == 1
    y0 = repmat(y0,1,W);
end
if size(n0,2) == 1
    n0 = repmat(n0,1,W);
end

%get matrices
M = array.M;
A = array.A;

%get phi
phi = zeros(array.Nn,size(x0,2));
for i = 1:L
   x = array.nodePosition(:,1);
   y = array.nodePosition(:,2);
   phi = phi + n0(i,:).*atan2(y-y0(i,:),x-x0(i,:));
end

%get th
areas = array.pathArea;
if size(f,1) == 1
    f = repmat(f,Np,1);
end
th = -M'*phi - 2*pi*A'*((A*A')\(areas.*f));

%refine th
if array.inductanceQ
    thp = round(th/2/pi)*2*pi;
    switch array.inductanceMode
        case 'uniform self'
            th = (th-thp)./(1+array.betaL.*array.Ic) + thp;
        case 'self'
            th = (th-thp)./(1+array.betaL.*array.Ic) + thp;
        case 'matrix'
            th = (spdiags(ones(Nj,1),0,Nj,Nj)+array.betaL*spdiags(array.Ic,0,Nj,Nj))\(th-thp) + thp;
    end
end

%get I
if array.customcpQ
    I = array.cp(array.Ic,th);
else
    I = array.Ic.*sin(th);
end

%get n
n = -A*round(th/2/pi);

%test if output n is equal to input n
nInput = zeros(Np,W);
for w = 1:W
    for i = 1:Np
        p = array.pathPosition{i};
        in = inpolygon(x0(:,w),y0(:,w),p(:,1),p(:,2));
        nInput(i,w) = nInput(i,w) + sum(n0(in,w));
    end
end
%if mean(mean(abs(n-nInput))) > 1E-10
%    warning('output n not equal to input n')
%end
end