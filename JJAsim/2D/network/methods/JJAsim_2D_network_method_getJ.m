function J = JJAsim_2D_network_method_getJ(array,I)
%J = JJAsim_2D_network_method_getJ(array,I)
%
%DESCRIPTION
% obtains the path current J from the junction current I.
%
%FIXED INPUT
% array     struct                  information about Josephson junction array    
% I         Nj* by W by Nt          junction current   
%
%OUTPUT
% J         Np by W by Nt           path current

if size(I,1) == 1
    I = repmat(I,array.Nj,1);
end
D = array.A*array.A';
sz = size(I);
J = D\(array.A*reshape(I,array.Nj,[]));
J = reshape(J,[array.Np,sz(2:end)]);
end