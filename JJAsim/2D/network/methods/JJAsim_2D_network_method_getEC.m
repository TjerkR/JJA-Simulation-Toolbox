function [EC,ECtot] = JJAsim_2D_network_method_getEC(array,V)
%[EC,ECtot] = JJAsim_2D_network_method_getEC(array,V)
%
%DESCRIPTION
% obtains the energy stored in the capacitors, called EC, from the 
% junction voltage V. EC lists it for all junctions individually, whereas 
% ECtot is the sum over all junctions.
%
%FIXED INPUT
% array     struct                  information about Josephson junction array    
% V         Nj by W by Nt           junction voltage
%
%OUTPUT
% EC        Nj by W by Nt           capacitance energy
% ECtot     W by Nt                 total capacitance energy of whole array

if array.capacitanceQ
    EC = 0.5*array.betaC.*V.^2;
else
    EC = zeros(size(V));
end
ECtot = shiftdim(sum(EC,1),1);
end