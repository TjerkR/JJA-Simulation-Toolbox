function [EM,EMtot] = JJAsim_2D_network_method_getEM(array,I)
%[EM,EMtot] = JJAsim_2D_network_method_getEM(array,I)
%
%DESCRIPTION
% obtains the energy stored in the magnetic self fields, called EM, from the 
% junction current I. EM lists it for all junctions individually, whereas 
% EMtot is the sum over all junctions.
%
%FIXED INPUT
% array     struct                  information about Josephson junction array    
% I         Nj by W by Nt           junction current
%
%OUTPUT
% EM        Nj by W by Nt           magnetic energy
% EMtot     W by Nt                 total magnetic energy of whole array
if size(I,1) ~= array.Nj
    error('length of first I dimension bust be Nj')
end
if array.inductanceQ
    sz = size(I);
    I = reshape(I,array.Nj,[]);
    switch array.inductanceMode
        case 'uniform self'
            EM = 0.5*array.betaL*I.^2;
        case 'self'
            EM = 0.5*array.betaL.*I.^2;
        case 'matrix'
            EM = I.*(array.betaL*I)/2;
    end
    EM = reshape(EM,[array.Nj,sz(2:end)]);
else
    EM = zeros(size(I));
end
EMtot = shiftdim(sum(EM,1),1);
end