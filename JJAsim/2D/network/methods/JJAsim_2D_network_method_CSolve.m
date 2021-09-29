function y = JJAsim_2D_network_method_CSolve(x,CComponentsReduced,islandComponents,nrOfComponents)
%y = JJAsim_2D_network_method_CSolve(x,CComponentsReduced,islandComponents,nrOfComponents)
%
%DESCRIPTION
% solves the system C*y = x, where C = M*M'. Takes into account that the network can have several 
% disconnectedcomponents and solves each component separately.
%
%FIXED INPUT
% x                     Nis by W                    right hand side of system
% CComponentsReduced    1 by nrOfComponents cell    C matrix for each network components
% islandComponents      Nis by 1                    enumerates to which component each island belongs
% nrOfComponents        1 by 1                      number of network components            
%
%OUTPUT
% y                     Nis by W                    solution obeying C*y = x.

W = size(x,2);
Z = zeros(1,W);
y = x;
for c = 1:nrOfComponents
    ind = islandComponents == c;
    xp = x(ind,:);
    xp = xp(1:end-1,:);
    y(ind,:) = cat(1,CComponentsReduced{c}\xp,Z);
end

