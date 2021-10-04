function array = generate_array(N, L, NameValueArgs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments
    N (1,1) double
    L (1,1) double
    NameValueArgs.hole_shape = 'square'
    NameValueArgs.Ic = 1
end

hole_shape = NameValueArgs.hole_shape;
Ic = NameValueArgs.Ic;


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

