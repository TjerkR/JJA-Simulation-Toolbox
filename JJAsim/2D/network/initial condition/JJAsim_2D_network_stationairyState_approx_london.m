function [th,I] = JJAsim_2D_network_stationairyState_approx_london(array,n,f)
%[th,I] = JJAsim_2D_network_stationairyState_approx_london(array,n,f)
%
%DESCRIPTION
% - Generate an approximate stationairy state for a 2D network with the 
%   within the London's gauge, which can be used as an initial guess for 
%   JJAsim_2D_network_stationairyState or initial condition for 
%   JJAsim_2D_network_simulate. 
% - The output is found the phase zone z = n. On can use
%   JJAsim_2D_network_method_changePhaseZone to convert it to other phase 
%   zones. 
%
%FIXED INPUT
% array      struct      information about Josephson junction array
% n          Np* by W    desired vortex configuration
% f          np* by W    frustration factor 
%
%OUTPUT
% th         Nj by W     approximate gauge invariant phase difference
% I          Nj by W     approximate junction current

if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

if mean(mean(abs(round(n)-n))) > 1E-10
    error('n must be integers');
end

%get sizes
Np = array.Np;
W = max([size(n,2),size(f,2)]);

%check input
n = JJAsim_method_checkInput(n,'double',[Np,W],[1,1],'n');
f = JJAsim_method_checkInput(f,'double',[Np,W],[1,1],'f');

%get matrices
A = array.A;

%get th
areas = array.pathArea;
if size(f,1) == 1
    f = repmat(f,Np,1);
end
df = 2*pi*(n - areas.*f);
th = A'*((A*A')\df);

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

%check n
if mean(mean(abs(A*round(th/2/pi)))) > 1E-10
    warning('the found output th does not have the specified vortex configuration')
end
end