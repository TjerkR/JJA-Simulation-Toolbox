function out = JJAsim_2D_network_stationairyState(array,inputMode,IExt,f,z,th1,varargin)
%out = JJAsim_2D_network_stationairyState(array,inputMode,IExt,f,z,th1,varargin)
%
%DESCRIPTION
% - Finds stationairy states for 2D josephson junction networks.
% - The array can be an arbitrary 2D network. This is specified by array, which is 
%   created with JJAsim_2D_network_create.
% - Requires an initial guess th1. This can be generated with the function call
%   JJAsim_2D_network_stationairyState_approx_london
% - A solution need not exist, or if it exits, it may not be found. solutionQ 
%   lists if a solution is found.
% - To find a solution, Newtons iteration is used. The number of iterations is 
%   stored in nrOfSteps. The maximum number of steps and acceptance tolerance 
%   can be specified with newton_steps and newton_tol respectively.
% - Several problems can be computed in parallel in one function call. inputMode
%   determines the way to specify these; the options are 'sweep' and 'list'.
% - All physical quantities are normalized, see documentation.
%
%FIXED INPUT
% array            	 struct       	Geometric information of the network. Create 
%                                   with JJAsim_2D_network_create.
% inputMode      	 string       	'sweep' or 'list'. Determines how to specify
%                                   multiple simultaneous problems. Effects IExt,f,z,n. 
%                                   For 'list' the number of problems is W, whereas
%                                 	for 'sweep', the number of problems W = WI*Wf.                       
% IExt (sweep)     	 WI by 1        External current
%      (list)    	 W* by 1      
% f (sweep)          Np* by Wf*     Frustration factor
%   (list)           Np* by W*      
% z (sweep)          Np* by Wf*     Phase zone
%   (list)           Np* by W*      
% th1                Nj* by W*      Initial guess for gauge invariant phase differences.    
%
%VARIABLE INPUT
% repetitions        1 by 1         Number of repetitions for the whole problem set
% newton_steps     	 1 by 1        	number of steps for newtons algorithm.
% newton_tol      	 1 by 1       	acceptance tolerance for newtons algorithm.
% storethQ       	 1 by 1       	Determines if th is stored in output.     
% storenQ            1 by 1         Determines if n is stored in output. 
% storeIQ            1 by 1         Determines if I is stored in output. 
% storeEQ            1 by 1         Determines if E is stored in output. 
% computePartitions	 1 by 1        	Splits the problem set into computePartitions which are 
%                                   computed consecutively. Use to reduce memory usage or
%                                   for parallel computation. 
% parallelQ      	 1 by 1       	If true, the partitions are computed in parallel on 
%                                   multiple cores.  
% cores          	 1 by 1       	Number of cores used by parallel computing. If 0, it 
%                                   is set to the amount of available cores.            
%
%OUTPUT
% out.solutionQ    	 W by 1      	true if valid solution is found
% out.nrOfSteps      W by 1         number of newton iterations done to find solution
% out.th         	 Nj by W    	gauge invariant phase difference of every junction
% out.n            	 Np by W      	vortex configuration
% out.I            	 Nj by W      	junction current
% out.E           	 W by 1       	total energy stored in stationairy state

%get variable input
inputParameters = {
    'repetitions'            1    
    'storethQ'               true
    'storenQ'                true
    'storeIQ'                true
    'storeEQ'                true
    'newton_steps'         	 20
    'newton_tol'           	 1E-7
    'computePartitions'      1
    'parallelQ'              false
    'cores'                  0
    };

options = JJAsim_method_parseOptions(inputParameters,varargin,'JJAsim_2D_network_stationairyState');
repetitions = options.repetitions;
storethQ = logical(options.storethQ);
storenQ = logical(options.storenQ);
storeIQ = logical(options.storeIQ);
storeEQ = logical(options.storeEQ);
newton_steps = options.newton_steps;
newton_tol = options.newton_tol;
computePartitions = options.computePartitions;
parallelQ = options.parallelQ;
cores = options.cores;

%get sizes
Nn = array.Nn;
Nj = array.Nj;
Np = array.Np;
M = array.M;
A = array.A;

%check if the array struct corresponds to that of a 2D network with free boundaries.
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not 2D network'); 
end

%determine physical quantities
IExtBaseJ = array.IExtBaseJ;
pathArea = array.pathArea;  
areaCompact = array.pathAreaCompactQ;
Ic = array.Ic; 
IcCompact = array.IcCompactQ;
cpQ = array.customcpQ;
if cpQ
    cp = array.cp; dcp = array.dcp; icp = array.icp; 
else
    cp = [];       dcp = [];        icp = [];       
end
if array.inductanceQ
    betaL = array.betaL;
    switch array.inductanceMode
        case 'uniform self'; LMode = 0;
        case 'self'; LMode = 1; betaL = spdiags(betaL,0,Nj,Nj);
        otherwise; LMode = 1;
    end
end

%determine number of parallel problems Wtot
switch inputMode
    case 'sweep'       
        probFlag = 2; 
        W = prod([size(IExt,1),max([size(f,2),size(z,2)])]); 
    case 'list'        
        probFlag = 1; 
        W = max([size(IExt,1),size(f,2),size(z,2),size(th1,2)]);
    otherwise
        error('input mode does not exist')
end
Wtot = W*repetitions;

%check input size and type
IExt = JJAsim_method_checkInput(IExt,'double',W,[probFlag,1],'IExt');
[f,fCompact] = JJAsim_method_checkInput(f,'double',[Np,W],[1,probFlag],'f');
[z,zCompact] = JJAsim_method_checkInput(z,'double',[Np,W],[1,probFlag],'z');
if sum(sum(z - round(z) ~= 0)) > 1
   error('z must be integers') 
end
if strcmp(inputMode,'sweep') && size(f,2) > 1 && size(z,2) > 1 && size(f,2) ~= size(z,2)
   error('size of z inconsistent with size of f'); 
end


%check initial condition input
th1 = JJAsim_method_checkInput(th1,'double',[Nj,W],[1,1],'th1');
if size(th1,1) == 1
    th1 = repmat(th1,Nj,1);
end

%check other input
repetitions = JJAsim_method_checkInput(repetitions,'double',1,1,'repetitions');
storethQ = JJAsim_method_checkInput(storethQ,'logical',1,1,'storeThQ');
storenQ = JJAsim_method_checkInput(storenQ,'logical',1,1,'storenQ');
storeIQ = JJAsim_method_checkInput(storeIQ,'logical',1,1,'storeIQ');
storeEQ = JJAsim_method_checkInput(storeEQ,'logical',1,1,'storeEJtotQ');
computePartitions = JJAsim_method_checkInput(computePartitions,'double',1,1,'computePartitions');
parallelQ = JJAsim_method_checkInput(parallelQ,'logical',1,1,'parallelQ');
cores = JJAsim_method_checkInput(cores,'double',1,1,'cores');
if cores == 0
    cores = feature('numcores');
end

%get lists which entry into the input variables to use for each problem
switch inputMode
    case 'sweep'
        sweepW = [size(IExt,1),max([size(f,2),size(z,2)])];
        IExtProblemList = JJAsim_method_sweep(size(IExt,1),1,sweepW,repetitions);
        fProblemList = JJAsim_method_sweep(size(f,2),2,sweepW,repetitions);
        zProblemList = JJAsim_method_sweep(size(z,2),2,sweepW,repetitions);
    case 'list'
        IExtProblemList = JJAsim_method_list(size(IExt,1),W,repetitions);
        fProblemList = JJAsim_method_list(size(f,2),W,repetitions);
        zProblemList = JJAsim_method_list(size(z,2),W,repetitions); 
end

th1ProblemList = JJAsim_method_list(size(th1,2),W,repetitions); 

%partition the problem lists in computePartitions partitions.
[IExtProblemSelection,IExtTableSelection,WIExt,partitionSize] = ...
    JJAsim_method_partition(IExtProblemList,computePartitions);
[fProblemSelection,fTableSelection,Wf] = JJAsim_method_partition(...
    fProblemList,computePartitions);
[zProblemSelection,zTableSelection,Wz] = JJAsim_method_partition(...
    zProblemList,computePartitions);
[th1ProblemSelection,th1TableSelection,Wth1] = JJAsim_method_partition(...
    th1ProblemList,computePartitions);

%initialize output variables
if storethQ;   	thOut = zeros(Nj,partitionSize,computePartitions); end
if storenQ;   	nOut = zeros(Np,partitionSize,computePartitions); end
if storeIQ;   	IOut = zeros(Nj,partitionSize,computePartitions); end
if storeEQ;     EOut = zeros(partitionSize,computePartitions); end
solutionQ = zeros(partitionSize,computePartitions);
nrOfSteps = zeros(partitionSize,computePartitions);

%set parallel pool
if parallelQ
    %initialize the parallel pool
    JJAsim_method_poolinit(cores);
else
    %remove any present parallel pool to make up free space.
    JJAsim_method_poolinit(1);
end

%start computation of all partitions
if parallelQ
    
    %distribute the partitions into runs where each run contains as much partitions as
    %there are cores. This is done such that for each run separately the selections can
    %be prepared beforehand and no unnecessary input data is send to the cores. 
    runs = ceil(computePartitions/cores);
    for P = 1:runs
        
        %prepare input data selection for this run
        pInd = (P-1)*cores + 1 : min(P*cores,computePartitions);
        coreInd = pInd - (P-1)*cores;
        IExtP = JJAsim_method_tocell(IExt,IExtTableSelection(:,pInd),1,2);
        fP = JJAsim_method_tocell(f,fTableSelection(:,pInd),2,2);
        zP = JJAsim_method_tocell(z,zTableSelection(:,pInd),2,2);
        th1P = JJAsim_method_tocell(th1,th1TableSelection(:,pInd),2,2);
        IprobP = IExtProblemSelection(:,pInd);
        fprobP = fProblemSelection(:,pInd);
        zprobP = zProblemSelection(:,pInd);
        th1probP = th1ProblemSelection(:,pInd);

        %prepare temporary output variables    
        if storethQ;  	thOutTemp = zeros(size(thOut,1),size(thOut,2),cores); end
        if storenQ;  	nOutTemp = zeros(size(nOut,1),size(nOut,2),cores); end
        if storeIQ;  	IOutTemp = zeros(size(IOut,1),size(IOut,2),cores); end
        if storeEQ;     EOutTemp = zeros(size(EOut,1),cores); end
        solutionQTemp = zeros(size(solutionQ,1),cores);
        nrOfStepsTemp = zeros(size(nrOfSteps,1),cores);
        
        %execute run
        if ~array.inductanceQ 
            %solve initial condition for array without inductance
            parfor p = coreInd
                [thOutTemp(:,:,p),nOutTemp(:,:,p),IOutTemp(:,:,p),EOutTemp(:,p),solutionQTemp(:,p),...
                    nrOfStepsTemp(:,p)] = ...
                    JJAsim_2D_network_stationairyState_priv_noscr_CPU(Nn,Nj,Np,M,A,...
                    partitionSize,IExtBaseJ,pathArea,areaCompact,IExtP{p},IprobP(:,p),...
                    WIExt(p),fP{p},fCompact,fprobP(:,p),Wf(p),...
                    zP{p},zCompact,zprobP(:,p),Wz(p),th1P{p},th1probP(:,p),Wth1(p),...
                    Ic,IcCompact,cpQ,cp,dcp,icp,storethQ,storenQ,storeIQ,storeEQ,newton_steps,...
                    newton_tol);
            end
        else 
            %solve initial condition for array with inductance
            parfor p = coreInd
                [thOutTemp(:,:,p),nOutTemp(:,:,p),IOutTemp(:,:,p),EOutTemp(:,p),solutionQTemp(:,p),...
                    nrOfStepsTemp(:,p)] = ...
                    JJAsim_2D_network_stationairyState_priv_scr_CPU(Nn,Nj,Np,M,A,...
                    betaL,LMode,partitionSize,IExtBaseJ,pathArea,areaCompact,IExtP{p},...
                    IprobP(:,p),WIExt(p),fP{p},fCompact,fprobP(:,p),Wf(p),...
                    zP{p},zCompact,zprobP(:,p),Wz(p),...
                    th1P{p},th1probP(:,p),Wth1(p),Ic,IcCompact,cpQ,cp,dcp,icp,...
                    storethQ,storenQ,storeIQ,storeEQ,newton_steps,newton_tol);
            end
        end
        
        %store data from this run
        if storethQ;  	thOut(:,:,pInd) = thOutTemp(:,:,coreInd); end
        if storenQ;  	nOut(:,:,pInd) = nOutTemp(:,:,coreInd); end
        if storeIQ;  	IOut(:,:,pInd) = IOutTemp(:,:,coreInd); end
        if storeEQ;     EOut(:,pInd) = EOutTemp(:,coreInd); end
        solutionQ(:,pInd) = solutionQTemp(:,coreInd);
        nrOfSteps(:,pInd) = nrOfStepsTemp(:,coreInd);
    end
else
    %execute simulation for all partitions on single core
    for p = 1:computePartitions
        
        %get data selection for this partition
        Iprob = IExtProblemSelection(:,p);    Isel = IExtTableSelection(:,p);
        fprob = fProblemSelection(:,p);       fsel = fTableSelection(:,p);
        zprob = zProblemSelection(:,p);       zsel = zTableSelection(:,p);
        th1prob = th1ProblemSelection(:,p);   th1sel = th1TableSelection(:,p);

        %execute simulation for this partition
        if ~array.inductanceQ
            %solve initial condition for array without inductance
            [thp,np,Ip,Ep,solutionQ(:,p),nrOfSteps(:,p)] = ...
                JJAsim_2D_network_stationairyState_priv_noscr_CPU(Nn,Nj,Np,M,A,...
                partitionSize,IExtBaseJ,pathArea,areaCompact,IExt(Isel,:),Iprob,...
                WIExt(p),f(:,fsel),fCompact,fprob,Wf(p),...
                z(:,zsel),zCompact,zprob,Wz(p),th1(:,th1sel),th1prob,Wth1(p),...
                Ic,IcCompact,cpQ,cp,dcp,icp,storethQ,storenQ,storeIQ,storeEQ,newton_steps,...
                newton_tol);
        else
            %solve initial condition for array with inductance
            [thp,np,Ip,Ep,solutionQ(:,p),nrOfSteps(:,p)] = ...
                JJAsim_2D_network_stationairyState_priv_scr_CPU(Nn,Nj,Np,M,A,...
                betaL,LMode,partitionSize,IExtBaseJ,pathArea,areaCompact,IExt(Isel,:),...
                Iprob,WIExt(p),f(:,fsel),fCompact,fprob,Wf(p),...
                z(:,zsel),zCompact,zprob,Wz(p),th1(:,th1sel),...
                th1prob,Wth1(p),Ic,IcCompact,cpQ,cp,dcp,icp,storethQ,storenQ,storeIQ,storeEQ,...
                newton_steps,newton_tol);
        end
        
        %store data
        if storethQ; thOut(:,:,p) = thp; end
        if storenQ; nOut(:,:,p) = np; end
        if storeIQ; IOut(:,:,p) = Ip; end
        if storeEQ; EOut(:,p) = Ep; end 
    end
end

%store output variables in out 
solutionQ = reshape(solutionQ,[],1);
solutionQ = solutionQ(1:Wtot);
out.solutionQ = solutionQ;
nrOfSteps = reshape(nrOfSteps,[],1);
nrOfSteps = nrOfSteps(1:Wtot);
out.nrOfSteps = nrOfSteps;

if storethQ
    thOut = reshape(thOut,Nj,partitionSize*computePartitions);
    out.th = thOut(:,1:Wtot); 
end
if storenQ
    nOut = reshape(nOut,Np,partitionSize*computePartitions);
    out.n = nOut(:,1:Wtot); 
end
if storeIQ
    IOut = reshape(IOut,Nj,partitionSize*computePartitions);
    out.I = IOut(:,1:Wtot); 
end
if storeEQ
    EOut = reshape(EOut,partitionSize*computePartitions,1);
    out.E = EOut(1:Wtot); 
end
end
