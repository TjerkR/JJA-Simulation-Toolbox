function array = JJAsim_2D_network_adjust(array,varargin)
%array = JJAsim_2D_network_adjust(array,varargin)
%
%DESCRIPTION
% Adjust physical properties of 2D electrical network of josephson junctions defined  
% in array. Cannot add/remove nodes or junctions.
%
%FIXED INPUT
% array               struct            josephson junction array
%
%OPTIONAL INPUT
% nodePosition   	  Nn by 2          (x,y) coordinates of the nodes
% IExtBase         	  Nn by 1           relative external current injected per node
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
%
%OUTPUT
% array               struct            adjusted array

%check if the array struct corresponds to that of a 2D network with free boundaries.
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

%parse input
inputParameters = {
    'nodePosition'        []
    'IExtBase'              []
    'Ic'                    []
    'customcpQ'             []
    'cp'                    []
    'dcp'                   []
    'icp'                   []
    'Ncp'                   []
    'Rn'                    []
    'betaC'                 []
    'betaL'                 []
    };
options = JJAsim_method_parseOptions(inputParameters,varargin,mfilename);
nodePosition = options.nodePosition;
IExtBase = options.IExtBase;
customcpQ = options.customcpQ;
cp = options.cp;
dcp = options.dcp;
icp = options.icp;
Ncp = options.Ncp;
Ic = options.Ic;
Rn = options.Rn;
betaC = options.betaC;
betaL = options.betaL;

Nn = array.Nn;
Nj = array.Nj;
Np = array.Np;
Nc = array.nrOfConnectedComponents;

%check IExtBase input
if ~isempty(IExtBase)
    
    %check IExtBase
    array.IExtBase = JJAsim_method_checkInput(IExtBase,'double',Nn,0,'IExtBase');
    for c = 1:Nc
        if sum(array.IExtBase(array.nodeComponents == c)) > 1E-10
            error('IExtBase must add up to 0 for each network component');
        end
    end
    
    %adjust junction IExt base
    Cred = array.CComponentsReduced;
    iscomp = array.nodeComponents;
    array.IExtBaseJ = array.M'*JJAsim_2D_network_method_CSolve(IExtBase,Cred,iscomp,Nc);
end

%check current-phase relation input
if customcpQ == false
    if ~isempty(cp)
        error('cp is specified while customcpQ is set to false')
    end
    if ~isempty(dcp)
        error('dcp is specified while customcpQ is set to false')
    end
    if ~isempty(icp)
        error('icp is specified while customcpQ is set to false')
    end
    
    if ~isempty(Ncp)
        error('Ncp is specified while customcpQ is set to false')
    end
end
if ~isempty(cp) || ~isempty(dcp) || ~isempty(icp) || ~isempty(Ncp)
    customcpQ = true;
end
if customcpQ 
    if isempty(customcpQ)
        error('customcpQ is not specified.')
    end
    if isempty(cp)
        error('cp is not specified.')
    end
    if isempty(dcp)
        error('dcp is not specified.')
    end
    if isempty(icp)
        error('icp is not specified.')
    end
    if isempty(Ncp)
        error('Ncp is not specified.')
    end
    array.customcpQ = true;
    array.cp = cp;
    array.dcp = dcp;
    array.icp = icp;
    array.Ncp = Ncp;
else
    Ncp = 1;
end

%check Ic input
if ~isempty(Ic) 
    array.Ic = Ic;
end
[array.Ic,array.IcCompactQ] = JJAsim_method_checkInput(array.Ic,'double',[Nj,Ncp],[1,0],'Ic');

%check Rn input
if ~isempty(Rn)
    [array.Rn,array.RnCompactQ] = JJAsim_method_checkInput(Rn,'double',Nj,1,'Rn');
    if sum(array.Rn < 1E-10) > 0
        error('normal state resistance of all junctions must be larger than zero');
    end
end

%check betaC input
if ~isempty(betaC)
    [array.betaC,array.betaCCompactQ] = JJAsim_method_checkInput(betaC,'double',Nj,1,'betaC');
    if sum(array.betaC < 0) > 0
        error('capacitance of all junctions must zero or positive');
    end
    if mean(abs(array.betaC)) < 1E-10
        array.capacitanceQ = false;
        array.betaC = 0;
    else
        array.capacitanceQ = true;
    end
end

%check nodePosition input
if ~isempty(nodePosition)
    array.nodePosition = JJAsim_method_checkInput(nodePosition,'double',[Nn,2],[0,0],'nodePosition');
    if size(unique(round(array.nodePosition,10),'rows'),1) < size(array.nodePosition,1)
        error('all nodes must have distinct coordinates');
    end
    
    %adjust junctionPosition
    juncis1x = array.nodePosition(array.junctionIsland1,1);
    juncis1y = array.nodePosition(array.junctionIsland1,2);
    juncis2x = array.nodePosition(array.junctionIsland2,1);
    juncis2y = array.nodePosition(array.junctionIsland2,2);
    array.junctionPosition = [juncis1x,juncis1y,juncis2x,juncis2y];
    
    %adjust path area, centre and node positions
    pathArea = zeros(Np,1);
    pathCentroid = zeros(Np,2);
    pathPosition = cell(Np,1);
    pathNodes = array.pathNodes;
    pathJunctions = array.pathJunctions;
    pathDirection = array.pathDirection;
    for np = 1:Np
        %load path
        pNodes = pathNodes{np};
        pJunc = pathJunctions{np};
        pDir = pathDirection{np};
        
        %calculate area and centroid
        [parea,pcx,pcy] = JJAsim_2D_method_getpolygon(array.nodePosition(pNodes,1),array.nodePosition(pNodes,2));
        
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
        pathPosition{np} = array.nodePosition(pathNodes{np},:);
    end
    array.pathNodes = pathNodes;
    array.pathJunctions = pathJunctions;
    array.pathDirection = pathDirection;
    array.pathArea = pathArea;
    array.pathAreaCompactQ = false;
    array.pathCentroid = pathCentroid;
    array.pathPosition = pathPosition;
    
    %reconstruct cycle space A, can be different because paths might be reversed.
    cycleBaseEntries = full(sum(sum(array.A ~= 0)));
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
end

%check betaL input
if ~isempty(betaL)
    if mean(abs(reshape(betaL,[],1))) < 1E-10
        array.inductanceQ = false;
        betaL = 0;
    else
        array.inductanceQ = true;
    end
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
end