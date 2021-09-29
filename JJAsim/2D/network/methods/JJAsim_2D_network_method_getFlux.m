function Phi = JJAsim_2D_network_method_getFlux(array,I)
%Phi = JJAsim_2D_network_method_getFlux(array,I)
%
%DESCRIPTION
% obtains the induced magnetic flux Phi from junction current I.
%
%INPUT
% array     struct                  information about Josephson junction array    
% I         Nj* by W by Nt          junction current   
%
%OUTPUT
% Phi       Np by W by Nt           induced magnetic flux

if size(I,1) == 1
    I = repmat(I,array.Nj,1);
end
sz = size(I);
I = reshape(I,array.Nj,[]);
if array.inductanceQ
    switch array.inductanceMode
        case 'uniform self'
            I = I*array.betaL;
        case 'self'
            I = I.*array.betaL;
        case 'matrix'
            I = array.betaL*I;
    end
else
    I = zeros([array.Nj,sz(2:end)]);
end
Phi = reshape(array.A*I,[array.Np,sz(2:end)])/2/pi;
end