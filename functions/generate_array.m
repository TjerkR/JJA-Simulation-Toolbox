function array = generate_array(N, L, hole_shape, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Ic = 1;
try
    if varargin{1} == 'Ic'
        Ic = varargin{2};
    end
end

%%%%% TEMP %%%%%%%%%%%%%%%%%%%%%%%%%
% Nj = 2*N*(N-1);
% variation = 0.1 * rand(Nj,1);
% Ic = 100*ones(Nj,1) + variation;

%%%%% TEMP %%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize array
array = JJAsim_2D_network_square(N,N,1,1,'y','Ic',Ic); %TODO: look into syntax
xy = array.nodePosition;

switch hole_shape
    case 'square'
        % L by L square
        bound1 = (N-1)/2 - L/2;
        bound2 = (N-1)/2 + L/2;
        bound3 = (N-1)/2 - L/2;
        bound4 = (N-1)/2 + L/2;

        nodeNrs = (xy(:,1) > bound1 & xy(:,1) < bound2) ...
            & (xy(:,2) > bound3 & xy(:,2) < bound4);
        nodeNrs = find(nodeNrs);
    case 'circle'
        % Circle with diameter L
        nodeNrs = (xy(:,1)-(N-1)/2).^2 + (xy(:,2)-(N-1)/2).^2 < (L/2)^2;
    case 'diamond'
        % Diamond with diagonal L parallel to array axes
        a = (N-1)/2;
        nodeNrs = abs(xy(:,1)-a) + abs(xy(:,2)-a) < L/2;
    otherwise
        error('Unknown hole shape')
end

% Remove nodes corresponding to hole
array = JJAsim_2D_network_removeNodes(array, nodeNrs);

end

