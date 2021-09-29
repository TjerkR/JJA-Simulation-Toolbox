function array = JJAsim_2D_network_square(Nx,Ny,ax,ay,currentDirection,varargin)
%array = JJAsim_2D_network_square(Nx,Ny,ax,ay,direction,varargin)
%
%DESCRIPTION
% generates a Nx by Ny square array. Inherits optional input from 
% JJAsim_2D_network_create.
%
%FIXED INPUT
% Nx                1 by 1   number of nodes in x-direction
% Ny                1 by 1   number of nodes in y-direction
% ax                1 by 1   node distance in x-direction
% ay                1 by 1   node distance in y-direction
% currentDirection  string   'x' or 'y'. Sets a homogeneous external 
%                            current in the respective  dimension.
%                           
%VARIABLE INPUT
% (inherited from JJAsim_2D_network_create)
%
%OUTPUT
% array            struct    container for all information of the array. 

if nargin < 2
    Ny = Nx;
end
if nargin < 3
    ax = 1;
end
if nargin < 4
    ay = 1;
end
if nargin < 5
    currentDirection = 'x';
end

xnr = repmat((1:Nx)',1,Ny);
ynr = repmat((1:Ny),Nx,1);
x = (xnr-1)*ax;
y = (ynr-1)*ay;
nodePosition = [reshape(x,[],1),reshape(y,[],1)];

jc = [1,1,2,1;1,1,1,2];
jcs = zeros(2*Nx*Ny,4);
for iy = 0:Ny-1
    for ix = 0:Nx-1
        ind = (1:2) + 2*(ix + iy*Nx);
        jcs(ind,:) = jc + [ix,iy,ix,iy;ix,iy,ix,iy];
    end
end

invalid = jcs(:,1) > Nx | jcs(:,1) < 1 | jcs(:,2) > Ny | jcs(:,2) < 1 |jcs(:,3) > Nx | jcs(:,3) < 1 | jcs(:,4) > Ny | jcs(:,4) < 1;
jcs(invalid,:) = [];

junctionIsland1 = jcs(:,1) + (jcs(:,2)-1)*Nx;
junctionIsland2 = jcs(:,3) + (jcs(:,4)-1)*Nx;

IExtBase = zeros(Nx*Ny,1);
switch currentDirection
    case 'x'
        IExtBase(xnr == 1) = 1;
        IExtBase(xnr == Nx) = -1;
    case 'y'
        IExtBase(ynr == 1) = 1;
        IExtBase(ynr == Ny) = -1;
    otherwise
        error('unrecognized direction');
end
array = JJAsim_2D_network_create(nodePosition,[junctionIsland1,junctionIsland2],IExtBase,varargin{:});
end