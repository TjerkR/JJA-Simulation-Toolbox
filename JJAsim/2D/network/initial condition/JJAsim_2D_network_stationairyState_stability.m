function [stableQ,eigv] = JJAsim_2D_network_stationaryState_stability(array,th,parallelQ)
%[stableQ,eigv] = JJAsim_2D_network_stationaryState_stability(array,th)
%
%DESCRIPTION
% - returns if a stationary state th is dynamically stable, defined
%   by the criterion that the maximum eigenvalue of the Jacobian of the
%   system is non-positive.
% - Can compute W problems in one function call.
% - Also outputs the eigenvalues of the Jacobian.
% - Does not explicitly check if th is a valid stationary state.
%
%FIXED INPUT
% array     struct      array defining a network of Josephson junctions.
% th        Nj by W     stationary state to check the stability of. 
%
%OUTPUT
% stableQ   1 by W      true if stationary state is stable
% eigv      Nj by W     listing all eigenvalues of the Jacobian matrix at th

if nargin < 3
    parallelQ = false;
end

%check if array is a 2D network
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

%get sizes
W = size(th,2);
Nj = array.Nj;

%check input
th = JJAsim_method_checkInput(th,'double',[Nj,W],[0,0],'th');

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
    if parallelQ
        parfor w = 1:W
            J = -diag(dth(:,w))-X;
            eigv(:,w) = real(eig(J));
            stableQ(w) = max(eigv(:,w),[],1) < 0;
        end
    else
        for w = 1:W
            J = -diag(dth(:,w))-X;
            eigv(:,w) = real(eig(J));
            stableQ(w) = max(eigv(:,w),[],1) < 0;
        end
    end
else
    M = array.M;
    if parallelQ
        parfor w = 1:W
            S = M*spdiags(dth(:,w),0,Nj,Nj);
            J = -M'*JJAsim_2D_network_method_CSolve(full(S),array.CComponentsReduced,...
                array.nodeComponents,array.nrOfConnectedComponents);
            eigv(:,w) = real(eig(J));
            stableQ(w) = max(eigv(:,w),[],1) < 1E-8;
        end
    else
        for w = 1:W
            S = M*spdiags(dth(:,w),0,Nj,Nj);
            J = -M'*JJAsim_2D_network_method_CSolve(full(S),array.CComponentsReduced,...
                array.nodeComponents,array.nrOfConnectedComponents);
            eigv(:,w) = real(eig(J));
            stableQ(w) = max(eigv(:,w),[],1) < 1E-8;
        end
    end
end
end