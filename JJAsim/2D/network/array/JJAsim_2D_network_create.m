function array = JJAsim_2D_network_create(nodePosition,junctionIslands,IExtBase,varargin)
%array = JJAsim_2D_network_create(nodePosition,junctionIslands,IExtBase,varargin)
%
%DESCRIPTION
% - Generate an 2D electrical network of Josephson junctions with arbitrary network 
%   structure. The resulting struct array is used in the JJAsim package to do 
%   simulations.
% - Includes geometrical information and physical quantities like critical current, 
%   resistance, capacitance and inductance. 
% - All physical quantities are dimensionless, see documentation for normalizations. 
% - Note that for networks which have the structure of a lattice one can use 
%   JJAsim_lattice_2D instead, as this special type generally allows for significantly
%   sped up computation. 
% - Inductances are specified as inductances between junctions (not between paths). 
%   Realistic values can be computed with the INDUCTsim package. 
%
%FIXED INPUT
% nodePosition        Nn by 2           x and y coordinates of the nodes of the array
% junctionIslands	  Nj by 2           node numbers of endpoints of each junction
% IExtBase         	  Nn by 1           relative external current injected per node
%
%OPTIONAL INPUT
% customcpQ           1 by 1            true if custom current-phase relation cp
% cp                  function_handle   custom current phase relation @(Ic,th) cp(Ic,th), 
%                                       which must be 2pi periodic.
% dcp                 function_handle   th-derivative of current phase relation
% icp                 function_handle   th-integral of current phase relation
% Ncp                 1 by 1            number of current-phase relation parameters
% Ic                  Nj* by Ncp        critical current of junctions
% Rn              	  Nj* by 1          normal state resistance of junctions
% betaC               Nj* by 1          capacitance of junctions (0 for no capacitance)
% betaL (self)        Nj* by 1          junction self inductance.           
%       (self+mutual) Nj by Nj          junction inductance matrix. Can be sparse or full. 
% cycleAlgorithm      string            'greedy' or 'spanningTree'. The algorithm used
%                                       to find the cycle space.
%
%OUTPUT:
% array               struct              container for all information of the array. 

%array struct with fields:
%   boudary                   string          boundary type (set to 'free')
%   type                      string          array type (set to 'network')
%   ndims                     1 by 1          number of spatial dimensions (set to 2)
%   Nn                        1 by 1          number of nodes
%   Nj                        1 by 1          number of junctions
%   Np                        1 by 1          number of paths (circuit rank of network).
%   N                         1 by 1          sqrt(Np)
%   nodePosition              Nn by 2         x and y coordinates of the nodes of the array
%   junctionIsland1           Nj by 1         node numbers of the junctions first node    
%   junctionIsland2           Nj by 1         node numbers of the junctions second node
%   junctionPosition          Nj by 4         Coordinates of endpoints of junctions, [x1,y1,x2,y2].
%   pathArea                  Np by 1         signed area of each path. Paths have positive 
%                                             orientation with right hand rule by convention.
%   pathCentroid              Np by 2         centroid of each path
%   pathPosition              Np by 1 cell    coordinates of the nodes in each path. each 
%                                             element is an Np by 2 array.
%   networkGraph              matlab graph    mathematical graph representation of the network
%   nrOfConnectedComponents   1 by 1          Abbrieviated Nc, number of connected components.          
%   M                         Nn by Nj        Kirchhoffs current rules
%   Mreduced                  Nn-1 by Nj      M with last row removed, the rows are linearly independent  
%   A                         Np by Nj        cycle (row) space or Kirchhoff's voltage rules.
%   nodeComponents            Nn by 1         (if Nc>1) to which component each node belongs    
%   junctionComponents        Nj by 1         (if Nc>1) to which component each junction belongs 
%   pathComponents            Np by 1         (if Nc>1) to which component each path belongs
%   networkGraphComponents    Nc by 1 cell    (if Nc>1) graph of each component   
%   MComponents               Nc by 1 cell    (if Nc>1) M matrix of each component
%   AComponents               Nc by 1 cell    (if Nc>1) A matrix of each component
%   IExtBase                  Nn by 1         relative external current injected per node
%   IExtBaseJ                 Nj by 1         MT*inv(M*MT)*IExtBase
%   Ic                        (Nj or 1) by 1  critical current of junctions
%   Rn                        (Nj or 1) by 1  normal state resistance of junctions
%   C                         (Nj or 1) by 1  capacitance of junctions (0 for no capacitance)
%   inductanceMode            string          'junction' of 'path'
%   L (junction)              Nj by 1         junction dependent self inductance.
%   L (path)                  Np by 1         path dependent self inductance. 
%                             or Np by Np     path inductance matrix. Can be sparse or full.

%import optional property-value pairs
inputParameters = {
    'Ic'                    1
    'customcpQ'             false
    'cp'                    @(Ic,th) Ic.*sin(th)
    'dcp'                   @(Ic,th) Ic.*cos(th)
    'icp'                   @(Ic,th) Ic.*(1-cos(th))
    'Ncp'                   1 
    'Rn'                    1
    'betaC'                 0
    'betaL'                 0
    'cycleAlgorithm'        'greedy'
    };
options = JJAsim_method_parseOptions(inputParameters,varargin,mfilename);
customcpQ = options.customcpQ;
cp = options.cp;
dcp = options.dcp;
icp = options.icp;
Ncp = options.Ncp;
Ic = options.Ic;
Rn = options.Rn;
betaC = options.betaC;
betaL = options.betaL;
cycleAlgorithm = options.cycleAlgorithm;

%array type
array.boundary = 'free';
array.type = 'network';
array.ndims = 2;

%get array sizes
Nn = size(nodePosition,1);
junctionIsland1 = junctionIslands(:,1);
junctionIsland2 = junctionIslands(:,2);
Nj = length(junctionIsland1);

%filter nodes that are disconnected
filtis = setdiff(1:Nn,unique([junctionIsland1;junctionIsland2]));
% %if ~isempty(filtis)
%     is = true(Nn,1)
%     is(filtis) = false
%     cis = cumsum(is)
%     junctionIsland1 = cis(junctionIsland1);
%     junctionIsland2 = cis(junctionIsland2);
%     Nn = cis(end);
%     nodePosition = nodePosition(is,:);
%     IExtBase = IExtBase(is);
%     warning('some nodes are disconnected, they are removed.') 
% end

%set array sizes
array.Nn = Nn;
array.Nj = Nj;

%check input
nodePosition = JJAsim_method_checkInput(nodePosition,'double',[Nn,2],[0,0],'nodePosition');
if size(unique(round(nodePosition,10),'rows'),1) < size(nodePosition,1)
    error('all nodes must have distinct coordinates');
end
junctionIsland1 = JJAsim_method_checkInput(junctionIsland1,'matlabInt',Nj,0,'junctionIsland1');
junctionIsland2 = JJAsim_method_checkInput(junctionIsland2,'matlabInt',Nj,0,'junctionIsland2');
array.IExtBase = JJAsim_method_checkInput(IExtBase,'double',Nn,0,'IExtBase');
array.customcpQ = customcpQ;
if customcpQ
    %build in check to see if cp and dcp are valid. 
    array.cp = cp;
    array.dcp = dcp;
    array.icp = icp;
    array.Ncp = Ncp;
else
    Ncp = 1;
end
[array.Ic,array.IcCompactQ] = JJAsim_method_checkInput(Ic,'double',[Nj,Ncp],[1,0],'Ic');
[array.Rn,array.RnCompactQ] = JJAsim_method_checkInput(Rn,'double',Nj,1,'Rn');
if sum(array.Rn < 1E-10) > 0
   error('normal state resistance of all junctions must be larger than zero'); 
end
[array.betaC,array.betaCCompactQ] = JJAsim_method_checkInput(betaC,'double',Nj,1,'betaC');
if sum(array.betaC < 0) > 0
   error('capacitance of all junctions must zero or positive'); 
end

if mean(abs(betaC)) < 1E-10
    array.capacitanceQ = false;
    array.betaC = 0;
else
    array.capacitanceQ = true;
end
if mean(abs(reshape(betaL,[],1))) < 1E-10
    array.inductanceQ = false;
    betaL = 0;
else
    array.inductanceQ = true;
end


%get kirchhoffs current rules
Mijk = zeros(2*Nj,3);
for i = 1:Nj
	Mijk(i,:) = [junctionIsland1(i),i,1];
	Mijk(i+Nj,:) = [junctionIsland2(i),i,-1];
end
M = sparse(Mijk(:,1),Mijk(:,2),Mijk(:,3),Nn,Nj);
array.M = M;

%make sure junctions dont start and end at the same node
if sum(junctionIsland1 == junctionIsland2) > 0
   error('junctions cannot start and end at the same node');
end

%get node and junction coordinates
array.nodePosition = nodePosition;
array.junctionIsland1 = junctionIsland1;
array.junctionIsland2 = junctionIsland2;
juncis1x = nodePosition(junctionIsland1,1);
juncis1y = nodePosition(junctionIsland1,2);
juncis2x = nodePosition(junctionIsland2,1);
juncis2y = nodePosition(junctionIsland2,2);
array.junctionPosition = [juncis1x,juncis1y,juncis2x,juncis2y];

%get graph of full network G
edgeTable = table([junctionIsland1,junctionIsland2],(1:Nj)','VariableNames',{'EndNodes','edgeNumbers'});
G = graph(edgeTable);
array.networkGraph = G;

%find connected components
nodeComponents = conncomp(G)';
nrOfConnectedComponents = max(nodeComponents);
array.nrOfConnectedComponents = nrOfConnectedComponents;

%get nr of paths (circuit rank)
Np = Nj - Nn + nrOfConnectedComponents;
array.Np = Np;
array.N = sqrt(Np);

%initialize component quantities
junctionComponents = zeros(Nj,1);
pathComponents = zeros(Np,1);
MComponents = cell(nrOfConnectedComponents,1);
AComponents = cell(nrOfConnectedComponents,1);
networkGraphComponents = cell(nrOfConnectedComponents,1);
size(nodeComponents);
%check if external current base adds up for each component
% for c = 1:nrOfConnectedComponents
%     sum(array.IExtBase(nodeComponents == c))
%     if sum(array.IExtBase(nodeComponents == c)) > 1E-10
%         error('IExtBase must add up to 0 for each network component');
%     end
% end

%initialize cycles
pathJunctions = cell(Np,1);
np = 1;
cycleBaseEntries = 0;

%compute all network components separately
for c = 1:nrOfConnectedComponents

    cIslands = nodeComponents  == c;
    cJunctions = logical(abs(M)'*cIslands ~= 0);
    junctionComponents(cJunctions) = c;
    
    %get network for this component
    Gcomp = subgraph(G,find(cIslands));
    networkGraphComponents{c} = Gcomp;
    
    %get MComponents
    MComponents{c} = M(cIslands,cJunctions);
    
    %find cycles
    switch cycleAlgorithm
        case 'greedy'
            for j = 1:numedges(Gcomp)
                nodes = Gcomp.Edges.EndNodes(1,:);
                edgeNr = Gcomp.Edges.edgeNumbers(1);
                Gcomp = rmedge(Gcomp,1);
                path = shortestpath(Gcomp,nodes(1),nodes(2));
                if ~isempty(path)
                    pathEdges = findedge(Gcomp,path(1:end-1),path(2:end));
                    pathEdges = [edgeNr;Gcomp.Edges.edgeNumbers(pathEdges)];
                    pathComponents(np) = c;
                    pathJunctions{np} = pathEdges;
                    cycleBaseEntries = cycleBaseEntries + length(pathEdges);
                    np = np + 1;
                end
            end
            
        case 'spanningTree'
            Tcomp = minspantree(Gcomp);
            edgesNotInTcomp = setdiff(Gcomp.Edges,Tcomp.Edges);
            edgesNotInTcompNodes = edgesNotInTcomp.EndNodes;
            edgesNotInTcompNumbers = edgesNotInTcomp.edgeNumbers;
            for p = 1:size(edgesNotInTcompNodes,1)
                path = shortestpath(Tcomp,edgesNotInTcompNodes(p,1),edgesNotInTcompNodes(p,2));
                pathEdges = findedge(Tcomp,path(1:end-1),path(2:end));
                pathEdges = [edgesNotInTcompNumbers(p);Tcomp.Edges.edgeNumbers(pathEdges)];
                pathComponents(np) = c;
                pathJunctions{np} = pathEdges;
                cycleBaseEntries = cycleBaseEntries + length(pathEdges);
                np = np + 1;
            end
            
        otherwise
            error('unrecognized cycleAlgorithm');
    end
end

%get reduced quantities used to solve y = M*M'\x when there are multiple components. 
CComponentsReduced = cell(1,nrOfConnectedComponents);
for c = 1:nrOfConnectedComponents
    Mr = MComponents{c};
    Mr = Mr(1:end-1,:);
    CComponentsReduced{c} = Mr*Mr';
end
array.CComponentsReduced = CComponentsReduced;

%get junction IExt base
array.IExtBaseJ = M'*JJAsim_2D_network_method_CSolve(IExtBase,CComponentsReduced,nodeComponents,nrOfConnectedComponents);

%get all information about the paths
pathNodes = cell(Np,1);
pathDirection = cell(Np,1);
for np = 1:Np
    pJunc = pathJunctions{np};
    [pathNodes{np},pathDirection{np}] = JJAsim_2D_network_pathproperties([junctionIsland1(pJunc),junctionIsland2(pJunc)]);
end

%construct path area, centre and node positions
pathArea = zeros(Np,1);
pathCentroid = zeros(Np,2);
pathPosition = cell(Np,1);
for np = 1:Np
    %load path
    pNodes = pathNodes{np};
    pJunc = pathJunctions{np};
    pDir = pathDirection{np};
    
    %calculate area and centroid
    [parea,pcx,pcy] = JJAsim_2D_method_getpolygon(nodePosition(pNodes,1),nodePosition(pNodes,2));
    
    %store area and centroid
    pathArea(np) = abs(parea);
    pathCentroid(np,:) = [pcx,pcy];
    
    %if area is negative, reverse the direction of the path. 
    if parea < 0
        pathNodes{np} = flipud(pNodes);
        pathJunctions{np} = flipud(pJunc);
        pathDirection{np} = -flipud(pDir);
    end

    %store node coordinates of paths 
    pathPosition{np} = nodePosition(pathNodes{np},:);
end
array.pathNodes = pathNodes;
array.pathJunctions = pathJunctions;
array.pathDirection = pathDirection;
array.pathArea = pathArea;
array.pathAreaCompactQ = false;
array.pathCentroid = pathCentroid;
array.pathPosition = pathPosition;

%construct cycle space A
Aijk = zeros(cycleBaseEntries,3);
AijkEntry = 1;
for np = 1:Np
    pJunc = pathJunctions{np};
    pDir = pathDirection{np}; 
    for p = 1:length(pJunc)
        Aijk(AijkEntry,:) = [np,pJunc(p),pDir(p)];
        AijkEntry = AijkEntry + 1;
    end
end
A = sparse(Aijk(:,1),Aijk(:,2),Aijk(:,3),Np,Nj);
array.A = A;

%get AComponents
for c = 1:nrOfConnectedComponents
    cJunctions = junctionComponents == c;
    cPaths = pathComponents == c;
    AComponents{c} = A(cPaths,cJunctions);
end

%if the number of connected components exceeds 1, store all related quantities.
array.nodeComponents = nodeComponents;
if nrOfConnectedComponents > 1
    array.junctionComponents = junctionComponents;
    array.pathComponents = pathComponents;
    array.MComponents = MComponents;
    array.AComponents = AComponents;
    array.networkGraphComponents = networkGraphComponents;
end

%check inductance
if size(betaL,1) == 1 && size(betaL,2) > 1
    error('betaL not the right size.');
end
betaL = JJAsim_method_checkInput(betaL,'double',[Nj,Nj],[1,1],'betaL');
if Nj ~= 1
    if size(betaL,1) == 1 && size(betaL,2) == 1
        if betaL < 0
            error('betaL matrix negative definite');
        end
        if betaL == 0
            array.inductanceQ = false;
        else
            array.inductanceMode = 'uniform self';
        end
    end
    if size(betaL,1) == Nj && size(betaL,2) == 1
        if min(betaL) < 0
            error('betaL matrix negative definite');
        end
        array.inductanceMode = 'self';
    end
    if size(betaL,1) == Nj && size(betaL,2) == Nj
        [~,p] = chol(betaL);
        if p ~= 0
            error('betaL matrix negative definite');
        end
        if sprank(A*betaL*A') < Np
            error('rank of A*betaL*AT is smaller than Np so it is not invertible. Likely some paths have no self-inductance, while others do, which can be done in principle but is not supported.')
        end
        array.inductanceMode = 'matrix';
    end
else
    if betaL < 0
        error('betaL matrix negative definite');
    end
    if betaL == 0
        array.inductanceQ = false;
    else
        array.inductanceMode = 'uniform self';
    end
end
array.betaL = betaL;
end
