function [stableQ,eigv] = JJAsim_2D_network_stationairyState_stability(array,th)
%[stableQ,eigv] = JJAsim_2D_network_stationairyState_stability(array,th)
%
%DESCRIPTION
% - returns if a stationairy state th is dynamically stable, defined
%   by the criterion that the maximum eigenvalue of the Jacobian of the
%   system is non-positive.
% - Can compute W problems in one function call.
% - Also outputs the eigenvalues of the Jacobian.
% - Does not explicitly check if th is a valid stationairy state.
%
%FIXED INPUT
% array     struct      array defining a network of Josephson junctions.
% th        Nj by W     stationairy state to check the stability of. 
%
%OUTPUT
% stableQ   1 by W      true if stationairy state is stable
% eigv      Nj by W     listing all eigenvalues of the Jacobian matrix at th

%check if array is a 2D network
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

W = size(th,2);
M = array.M;

%compute gradient of th vector
if array.customcpQ
    dth = array.dcp(array.Ic,th);
else
    dth = array.Ic.*cos(th);
end

%compute Jacobian, eigenvalues and stability
eigv = zeros(array.Nj,W);
stableQ = zeros(1,W);
if array.inductanceQ
    A = array.A;
    B = A*(array.betaL.*A');
    X = A'*(B\A);
    for w = 1:W
        J = -diag(dth(:,w))-X;
        eigv(:,w) = real(eig(J));
        stableQ(w) = max(eigv(:,w),[],1) < 0;
    end
else
    for w = 1:W
        J = -M'*JJAsim_2D_network_method_CSolve(array.M*diag(dth(:,w)),array.CComponentsReduced,...
        array.nodeComponents,array.nrOfConnectedComponents);
        eigv(:,w) = real(eig(J));
        stableQ(w) = max(eigv(:,w),[],1) < 1E-8;
    end
end
end