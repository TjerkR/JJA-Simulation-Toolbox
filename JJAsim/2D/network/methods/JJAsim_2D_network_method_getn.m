function n = JJAsim_2D_network_method_getn(array,th,z)
%n = JJAsim_2D_network_method_getn(array,th,z)
%
%DESCRIPTION
% obtains the vortex configuration n from phase difference th and phase zone z.
%
%INPUT
% array     struct                  information about Josephson junction array
% th        Nj* by W* by Nt         gauge invariant phase difference     
% z         Np* by W*               phase zone
%
%OUTPUT
% n         Np by W by Nt           vortex configuration

if size(th,1) == 1
    th = repmat(th,array.Nj,1);
end
sz = size(th);
n = array.A*reshape(round(th/2/pi),array.Nj,[]);
n = z - reshape(n,[array.Np,sz(2:end)]);
end