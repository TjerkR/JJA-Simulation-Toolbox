function U = JJAsim_2D_network_method_getU(array,V)
%U = JJAsim_2D_network_method_getU(array,V)
%
%DESCRIPTION
% obtains the node potential U from the junction voltage drop V. 
%
%FIXED INPUT
% array     struct                  information about Josephson junction array    
% V         Nj* by W by Nt          junction voltage drop
%
%OUTPUT
% U         Nn by W by Nt           node potential

if size(V,1) == 1
    V = repmat(V,array.Nj,1);
end
Cred = array.CComponentsReduced;
Ncomp = array.nrOfConnectedComponents;
if Ncomp > 1
    iscomp = array.nodeComponents;
else
    iscomp = ones(array.Nn,1);
end
sz = size(V);
U = -JJAsim_2D_network_method_CSolve(array.M*reshape(V,array.Nj,[]),Cred,iscomp,Ncomp);
U = reshape(U,[array.Nn,sz(2:end)]);
end