function array = JJAsim_2D_network_removeIslands(array,islandNrs)
%array = JJAsim_2D_network_removeIslands(array,islandNrs)
%
%DESCRIPTION
% remove islands from an existing 2D network of josephson junctions.
%
%FIXED INPUT
% array         struct            josephson junction array
% islandNrs     N by 1            list of island numbers to remove.
%
%OUTPUT
% array         struct            adjusted array

%check if the array struct corresponds to that of a 2D network with free boundaries.
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

%remove islands and junctions
islands = true(array.Nis,1);
islands(islandNrs) = false;
cislands = cumsum(islands);
islandPosition = array.islandPosition(islands,:);
j1 = array.junctionIsland1;
j2 = array.junctionIsland2;
ind = islands(j1) & islands(j2);
junctionIsland1 = cislands(j1(ind));
junctionIsland2 = cislands(j2(ind));
IExtBase = array.IExtBase(islands);
Nis = sum(islands);
Nj = length(junctionIsland1);

%filter islands that are disconnected
filtis = setdiff(1:Nis,unique([junctionIsland1;junctionIsland2]));
if ~isempty(filtis)
    is = true(Nis,1);
    is(filtis) = false;
    cis = cumsum(is); 
    junctionIsland1 = cis(junctionIsland1);
    junctionIsland2 = cis(junctionIsland2);
    Nis = cis(end);
    islandPosition = islandPosition(is,:);
    IExtBase = IExtBase(is);
    warning('some islands are disconnected, they are removed.') 
end

%set array sizes
array.Nis = Nis;
array.Nj = Nj;

%contract physical parameters
if ~array.IcCompactQ(1)
   array.Ic = array.Ic(ind,:);
end
if ~array.RnCompactQ
    array.Rn = array.Rn(ind);
end
if ~array.betaCCompactQ
    array.betaC = array.betaC(ind);
end
if array.inductanceQ
    switch array.inductanceMode
        case 'self'
            array.betaL = array.betaL(ind);
        case 'matrix'
            array.betaL = array.betaL(ind,ind);
            [~,p] = chol(array.betaL);
            if p ~= 0
                error('betaL matrix negative definite');
            end
    end
end

%recalculate kirchhoffs current rules
Mijk = zeros(2*Nj,3);
for i = 1:Nj
	Mijk(i,:) = [junctionIsland1(i),i,1];
	Mijk(i+Nj,:) = [junctionIsland2(i),i,-1];
end
M = sparse(Mijk(:,1),Mijk(:,2),Mijk(:,3),Nis,Nj);
array.M = M;

%get island and junction coordinates
array.islandPosition = islandPosition;
array.junctionIsland1 = junctionIsland1;
array.junctionIsland2 = junctionIsland2;
array.junctionPosition = array.junctionPosition(ind,:);

%get graph of full network G
edgeTable = table([junctionIsland1,junctionIsland2],(1:Nj)','VariableNames',{'EndNodes','edgeNumbers'});
G = graph(edgeTable);
array.networkGraph = G;

%find connected components
islandComponents = conncomp(G)';
if Nis > 0
    nrOfConnectedComponents = max(islandComponents);
else
    nrOfConnectedComponents = 0;
end
array.nrOfConnectedComponents = nrOfConnectedComponents;

%get nr of paths (circuit rank)
Np = Nj - Nis + nrOfConnectedComponents;
array.Np = Np;
array.N = sqrt(Np);

%initialize component quantities
junctionComponents = zeros(Nj,1);
pathComponents = zeros(Np,1);
MComponents = cell(nrOfConnectedComponents,1);
AComponents = cell(nrOfConnectedComponents,1);
networkGraphComponents = cell(nrOfConnectedComponents,1);

%initialize cycles
pathJunctions = cell(Np,1);
np = 1;
cycleBaseEntries = 0;

%compute all network components separately
for c = 1:nrOfConnectedComponents

    cIslands = islandComponents  == c;
    cJunctions = logical(abs(M)'*cIslands ~= 0);
    junctionComponents(cJunctions) = c;
    
    %get network for this component
    Gcomp = subgraph(G,find(cIslands));
    networkGraphComponents{c} = Gcomp;
    
    %get MComponents
    MComponents{c} = M(cIslands,cJunctions);
    
    %find cycles
    for j = 1:numedges(Gcomp)
        islands = Gcomp.Edges.EndNodes(1,:);
        edgeNr = Gcomp.Edges.edgeNumbers(1);
        Gcomp = rmedge(Gcomp,1);
        path = shortestpath(Gcomp,islands(1),islands(2));
        if ~isempty(path)
            pathEdges = findedge(Gcomp,path(1:end-1),path(2:end));
            pathEdges = [edgeNr;Gcomp.Edges.edgeNumbers(pathEdges)];
            pathComponents(np) = c;
            pathJunctions{np} = pathEdges;
            cycleBaseEntries = cycleBaseEntries + length(pathEdges);
            np = np + 1;
        end
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

%get all information about the paths
pathNodes = cell(Np,1);
pathDirection = cell(Np,1);
for np = 1:Np
    pJunc = pathJunctions{np};
    [pathNodes{np},pathDirection{np}] = JJAsim_2D_network_pathproperties([junctionIsland1(pJunc),junctionIsland2(pJunc)]);
end

%construct path area, centre and island positions
pathArea = zeros(Np,1);
pathCentroid = zeros(Np,2);
pathPosition = cell(Np,1);
for np = 1:Np
    %load path
    pNodes = pathNodes{np};
    pJunc = pathJunctions{np};
    pDir = pathDirection{np};
    
    %calculate area and centroid
    [parea,pcx,pcy] = JJAsim_2D_method_getpolygon(islandPosition(pNodes,1),islandPosition(pNodes,2));
    
    %store area and centroid
    pathArea(np) = abs(parea);
    pathCentroid(np,:) = [pcx,pcy];
    
    %if area is negative, reverse the direction of the path. 
    if parea < 0
        pathNodes{np} = flipud(pNodes);
        pathJunctions{np} = flipud(pJunc);
        pathDirection{np} = -flipud(pDir);
    end

    %store island coordinates of paths 
    pathPosition{np} = islandPosition(pathNodes{np},:);
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
array.islandComponents = islandComponents;
if nrOfConnectedComponents > 1
    array.junctionComponents = junctionComponents;
    array.pathComponents = pathComponents;
    array.MComponents = MComponents;
    array.AComponents = AComponents;
    array.networkGraphComponents = networkGraphComponents;
end

%check if betaL is positive definit
if size(array.betaL,1) == Nj && size(array.betaL,2) == Nj
    if sprank(A*array.betaL*A') < Np
        error('rank of A*betaL*AT is smaller than Np so it is not invertible. Likely some paths have no self-inductance, while others do, which can be done in principle but is not supported.')
    end
end

%adjust external current base
pind = IExtBase > 0;
mind = IExtBase < 0;
rescaleFlag = false;
zeroFlag = false;
if array.nrOfConnectedComponents > 1
    for c = 1:array.nrOfConnectedComponents
        ind = array.islandComponents == c;
        if sum(ind & mind) > 0 && sum(ind & pind) > 0
            Iscale = (sum(IExtBase(ind & pind))- sum(IExtBase(ind & mind)))/2;
            rescaleFlag = true;
            IExtBase(ind & mind) = -Iscale*IExtBase(ind & mind)/sum(IExtBase(ind & mind));
            IExtBase(ind & pind) = Iscale*IExtBase(ind & pind)/sum(IExtBase(ind & pind));
        else
            zeroFlag = true;
            IExtBase(ind & mind) = 0;
            IExtBase(ind & pind) = 0;
        end
    end
else
    if sum(mind) > 0 && sum(pind) > 0
        Iscale = (sum(IExtBase(pind))-sum(IExtBase(mind)))/2;
        rescaleFlag = true;
        IExtBase(mind) = -Iscale*IExtBase(mind)/sum(IExtBase(mind));
        IExtBase(pind) = Iscale*IExtBase(pind)/sum(IExtBase(pind));
    else
        zeroFlag = true;
        IExtBase(mind) = 0;
        IExtBase(pind) = 0;
    end
end
if zeroFlag && ~rescaleFlag
   warning('Some IExtBase entries are set to zero') 
end
if rescaleFlag && ~zeroFlag
   warning('Some IExtBase entries are rescaled') 
end
if rescaleFlag && zeroFlag
   warning('Some IExtBase entries are rescaled and others are set to zero') 
end
array.IExtBase = IExtBase;

%get junction IExt base
array.IExtBaseJ = M'*JJAsim_2D_network_method_CSolve(IExtBase,CComponentsReduced,islandComponents,nrOfConnectedComponents);


end