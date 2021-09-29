function array = JJAsim_2D_network_honeycomb(Nx,Ny,ax,ay,currentDirection,varargin)
%array = JJAsim_2D_network_honeycomb(Nx,Ny,ax,ay,direction,varargin)
%
%DESCRIPTION
% generates a Nx by Ny unit cells honeycomb array. A unit cell contains
% 4 nodes. Inherits optional input from JJAsim_2D_network_create.
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

xnr = repmat((1:Nx)',1,Ny,4);
ynr = repmat((1:Ny),Nx,1,4);
incellnr = repmat(shiftdim(1:4,-1),Nx,Ny,1);
x = ((xnr-1)*3 + shiftdim([0,0.5,1.5,2],-1))*ax;
y = (ynr-1 + shiftdim([0,1,1,0]/2,-1))*ay*sqrt(3);
nodePosition = [reshape(permute(x,[3,1,2]),[],1),reshape(permute(y,[3,1,2]),[],1)];

jc = [
    1 1 1  1 1 2
    1 1 2  1 2 1
    1 1 2  1 1 3
    1 1 3  1 2 4
    1 1 3  1 1 4
    1 1 4  2 1 1
    ];
jcs = zeros(6*Nx*Ny,6);
for iy = 0:Ny-1
    for ix = 0:Nx-1
        ind = (1:6) + 6*(ix + iy*Nx);
        jcs(ind,:) = jc + [ix,iy,0,ix,iy,0];
    end
end

invalid = jcs(:,1) > Nx | jcs(:,1) < 1 | jcs(:,2) > Ny | jcs(:,2) < 1 |jcs(:,4) > Nx | jcs(:,4) < 1 | jcs(:,5) > Ny | jcs(:,5) < 1;
jcs(invalid,:) = [];

junctionIsland1 = jcs(:,3) + 4*(jcs(:,1)-1) + 4*Nx*(jcs(:,2)-1);
junctionIsland2 = jcs(:,6) + 4*(jcs(:,4)-1) + 4*Nx*(jcs(:,5)-1);

IExtBase = zeros(Nx,Ny,4);
switch currentDirection
    case 'x'
        IExtBase(xnr == 1 & incellnr == 1) = 1;
        IExtBase(xnr == Nx  & incellnr == 4) = -1;
    case 'y'
        IExtBase(ynr == 1  & (incellnr == 1 | incellnr == 4)) = 1;
        IExtBase(ynr == Ny & (incellnr == 2 | incellnr == 3)) = -1;
    otherwise
        error('unrecognized direction');
end
IExtBase = reshape(permute(IExtBase,[3,1,2]),[],1);
array = JJAsim_2D_network_create(nodePosition,[junctionIsland1,junctionIsland2],IExtBase,varargin{:});
end