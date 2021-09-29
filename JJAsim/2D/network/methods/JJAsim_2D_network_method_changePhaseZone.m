function th_new = JJAsim_2D_network_method_changePhaseZone(array,th_old,z_old,z_new)
%th_new = JJAsim_2D_network_method_changePhaseZone(array,th_old,z_old,z_new)
%
%DESCRIPTION
% Expresses th in a different phase zone. 
%
%FIXED INPUT
% array     struct              information about Josephson junction array
% th_old    Nj* by W* by Nt   	gauge invariant phase difference    
% z_old     Np* by W*           phase zone of th_old
% z_new     Np* by W*           new phase zone to express th_old in
%
%OUTPUT
% th_new    Nj by W by Nt       gauge invariant phase difference in new phase zone

dz = z_new-z_old;
if size(dz,1) == 1
    dz = repmat(dz,array.Np,1);
end
dth = array.A\dz;
if mean(mean(abs(round(dth)-dth),1)) > 1E-10
    error('failed to convert th_old to the desired phase zone z_new.')
end
th_new = th_old + 2*pi*dth;
end