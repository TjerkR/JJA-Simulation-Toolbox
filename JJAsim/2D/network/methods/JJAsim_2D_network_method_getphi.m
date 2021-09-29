function phi = JJAsim_2D_network_method_getphi(array,th)
%phi = JJAsim_2D_network_method_getphi(array,th)
%
%DESCRIPTION
% obtains the gauge dependent phase at each node from the gauge invariant 
% phase difference th.
%
%FIXED INPUT
% array     struct               information about Josephson junction array    
% th        Nj* by W by Nt       gauge invariant phase difference
%
%OUTPUT
% phi       Nn by W by Nt        gauge dependent phase at each node

if size(th,1) == 1
    th = repmat(th,array.Nj,1);
end
Cred = array.CComponentsReduced;
Ncomp = array.nrOfConnectedComponents;
iscomp = array.nodeComponents;
sz = size(th);
phi = -JJAsim_2D_network_method_CSolve(array.M*reshape(th,array.Nj,[]),Cred,iscomp,Ncomp);
phi = reshape(phi,[array.Nn,sz(2:end)]);
end