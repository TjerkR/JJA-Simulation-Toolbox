function [EJ,EJtot] = JJAsim_2D_network_method_getEJ(array,th)
%[EJ,EJtot] = JJAsim_2D_network_method_getEJ(array,th)
%
%DESCRIPTION
% obtains energy stored in the Josephson elements, called EJ, from the gauge 
% invariant phase difference th. EJ lists it for all junctions individually,
% whereas EJtot is the sum over all junctions.
%
%FIXED INPUT
% array     struct                  information about Josephson junction array    
% th        Nj by W by Nt           gauge invariant phase difference
%
%OUTPUT
% EJ        Nj by W by Nt           Josephson energy
% EJtot     W by Nt                 total Josephson energy of whole array
if array.customcpQ
    EJ = array.icp(array.Ic,th);
else
    EJ = array.Ic.*(1-cos(th));
end
EJtot = shiftdim(sum(EJ,1),1);
end