function out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1,varargin)
%out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1,varargin)
%
%DESCRIPTION
% - Compute time evolution on 2D josephson junction networks. The array can be
%   an arbitrary 2D network. This is specified by array, which is created with 
%   JJAsim\_2D\_network\_create.
% - Requires an initial condition th1. A reasonable initial condition can be 
%   generated with the JJAsim\_2D\_network\_stationairyState\_approx functions.
% - several problems can be computed in one function call. 
% - if the array contains capacitance (i.e. array.capacitanceQ is true), the 
%   second timestep must also be specified with th2. 
% - All physical quantities are dimensionless, see documentation.
%
%FIXED INPUT
% array            	struct           	 Geometric information of the network, 
%                                     	 Create with JJAsim_2D_network_create.       
% t                 1 by Nt          	 simulated time points. 
% inputMode       	string            	 'sweep' or 'list'. Determines how to
%                                     	 specify multiply simultaneous problems.
%                                        Effects IExt, f, z and T. For 'list', the
%                                        number of parallel problems is W. For
%                                     	 'sweep', W = WI*Wf*WT.
% IExt (sweep)      WI by Nt*            External current
%      (list)       W* by Nt*   
% T    (sweep)      WT by Nt*            temperature parameter
%      (list)       W* by Nt* 
% f    (sweep)      Np* by Wf*           frustration factor
%      (list)       Np* by W*   
% z    (sweep)      Np* by Wf*           phase zone
%      (list)       Np* by W*    
% th1               Nj* by W*            Initial condition phase differences  
%
%VARIABLE INPUT
% th2             	Nj* by W*            phase differences at second timestep. Only if 
%                                    	 array.capacitanceQ. Equal to th1 if empty. 
% Vscheme       	string           	 Scheme of the voltage derivative dth/dt. Can
%                                   	 be set to 'forward', 'central' or 'backward'.
%                                    	 Only used if array.capacitanceQ = true. 
% storethQ          1 by 1             	 If true, phase difference th is stored in output.
% storeIQ           1 by 1             	 If true, junction current I is stored in output.
% storeVQ           1 by 1             	 If true, junction voltage V is stored in output.
% storeVtotQ        1 by 1             	 If true, total voltage Vtot is stored in output.
% storeEQ           1 by 1             	 If true, total energy E is stored in output.
% thTimePoints  	Nt by 1           	 Lists which for which timepoints th is stored. 
% ITimePoints       Nt by 1           	 Lists which for which timepoints I is stored.   
% VTimePoints       Nt by 1           	 Lists which for which timepoints V is stored. 
% VtotTimePoints  	Nt by 1           	 Lists which for which timepoints Vtot is stored. 
% ETimePoints       Nt by 1           	 Lists which for which timepoints E is stored. 
% repetitions      	1 by 1            	 Number of repetition of the whole problem set.
%                                        rep for short. The total nr of problems is Wtot = W*rep.
% computePartitions	1 by 1          	 Splits the problems into computePartitions which 
%                                     	 are computed separately. Use to reduce memory or
%                                        for parallel computation.
% parallelQ         1 by 1            	 If true, the partitions are computed in parallel 
%                                   	 on multiple cores.  
% cores           	1 by 1          	 Number of cores used by parallel computing. If 0, 
%                                     	 it is set to the amount of available cores.      
%
%OUTPUT:
% out.th          	Nj by Wtot by Ntth	 gauge invariant phase difference output
% out.I          	Nj by Wtot by NtI 	 junction current output
% out.V          	Nj by Wtot by NtV 	 junction voltage drop output
% out.Vtot        	Wtot by NtVtot    	 total voltage drop output
% out.E          	Wtot by NtE       	 total energy output

%get sizes
Nn = array.Nn;
Nj = array.Nj;
Np = array.Np;
Nt = length(t);

%optional input parameters and their default values
inputParameters = {
    'th2'                   []
    'Vscheme'               'forward'
    'thTimePoints'          true(Nt,1)
    'ITimePoints'           true(Nt,1)
    'VTimePoints'           true(Nt,1)
    'VtotTimePoints'        true(Nt,1)
    'ETimePoints'           true(Nt,1) 
    'storethQ'              true  
    'storeIQ'               true
    'storeVQ'               true   
    'storeVtotQ'            true
    'storeEQ'               true
    'repetitions'           1
    'computePartitions'     1
    'parallelQ'             false
    'cores'                 0
    };

options = JJAsim_method_parseOptions(inputParameters,varargin,'JJAsim_2D_network_simulate');

th2 = options.th2;
Vscheme = options.Vscheme;
thTimePoints = options.thTimePoints;
ITimePoints = options.ITimePoints;
VTimePoints = options.VTimePoints;
VtotTimePoints = options.VtotTimePoints;
ETimePoints = options.ETimePoints;
storethQ = logical(options.storethQ);
storeIQ = logical(options.storeIQ);
storeVQ = logical(options.storeVQ);
storeVtotQ = logical(options.storeVtotQ);
storeEQ = logical(options.storeEQ);
repetitions = options.repetitions;
computePartitions = options.computePartitions;
parallelQ = options.parallelQ;
cores = options.cores;

%check if the array struct corresponds to that of a 2D network with free boundaries.
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not a 2D network'); 
end

%get array quantities
M = array.M;
A = array.A;
IExtBaseJ = array.IExtBaseJ;
IExtBaseTotal = sum(array.IExtBase(array.IExtBase > 0));
pathArea = array.pathArea;  areaCompact = array.pathAreaCompactQ;
if Np > 0
if pathArea(1) == 0
    error('area of first path cannot be zero as it is used for normalization');
end
end
Rn = array.Rn;  RnCompact = array.RnCompactQ;
Ic = array.Ic;  IcCompact = array.IcCompactQ;
cpQ = array.customcpQ;
if cpQ; cp =array.cp; icp = array.icp;  else; cp = []; icp = []; end
if array.inductanceQ 
    betaL = array.betaL;
end
if array.inductanceQ
    switch array.inductanceMode
        case 'uniform self'; LMode = 0;
        case 'self'; LMode = 1; betaL = spdiags(betaL,0,Nj,Nj);
        otherwise; LMode = 1;
    end
end

%th2 default value
if array.capacitanceQ
    if isempty(th2)
        th2 = th1;
    end
    switch Vscheme
        case 'forward'; Vscheme = 0;
        case 'central'; Vscheme = 1;
        case 'backward'; Vscheme = 2;
        otherwise; error('unrecognized Vscheme type')
    end
    betaC = array.betaC;
    betaCCompact = array.betaCCompactQ;
    
else
    th2 = [];
end

%check if time is linear if capacitance
if array.capacitanceQ
   if mean(abs(diff(t)-mean(diff(t)))) > 1E-8
       error('time must be linearly increasing if the array contains capacitance');
   end
end
dt = t(2)-t(1);

%determine W
switch inputMode
    case 'sweep';       probFlag = 2; W = prod([size(IExt,1),size(T,1),max(size(f,2),size(z,2))]);
    case 'list';        probFlag = 1; W = max([size(IExt,1),size(T,1),size(f,2),size(z,2),size(th1,2),size(th2,2)]);
    otherwise;          error('input mode does not exist')  
end
Wtot = W*repetitions;

%check input size and type
JJAsim_method_checkInput(t,'double',Nt,0);
[IExt,IExtCompact] = JJAsim_method_checkInput(IExt,'double',[W,Nt],[probFlag,1],'IExt');
[T,TCompact] = JJAsim_method_checkInput(T,'double',[W,Nt],[probFlag,1],'T');
[f,fCompact] = JJAsim_method_checkInput(f,'double',[Np,W],[1,probFlag],'f');
[z,zCompact] = JJAsim_method_checkInput(z,'double',[Np,W],[1,probFlag],'z');
if sum(sum(z - round(z) ~= 0)) > 1
   error('z must be integers') 
end
% if strcmp(inputMode,'sweep') && size(f,2) > 1 && size(f,2) ~= size(z,2)
%     error('size of z inconsistent with size of f');
% end

%check initial condition input
th1 = JJAsim_method_checkInput(th1,'double',[Nj,W],[1,1],'th1');
if size(th1,1) == 1
    th1 = repmat(th1,Nj,1);
end
if array.capacitanceQ
    th2 = JJAsim_method_checkInput(th2,'double',[Nj,W],[1,1],'th2_th');
    if size(th2,1) == 1
        th2 = repmat(th2,Nj,1);
    end
end

%check other input variables
thTimePoints = JJAsim_method_checkInput(thTimePoints,'logical',Nt,0,'thTimePoints');
ITimePoints = JJAsim_method_checkInput(ITimePoints,'logical',Nt,0,'ITimePoints');
VTimePoints = JJAsim_method_checkInput(VTimePoints,'logical',Nt,0,'VTimePoints');
VtotTimePoints = JJAsim_method_checkInput(VtotTimePoints,'logical',Nt,0,'VtotTimePoints');
ETimePoints = JJAsim_method_checkInput(ETimePoints,'logical',Nt,0,'ETimePoints');
storethQ = JJAsim_method_checkInput(storethQ,'logical',1,1,'storethQ');
storeIQ = JJAsim_method_checkInput(storeIQ,'logical',1,1,'storeIQ');
storeVQ = JJAsim_method_checkInput(storeVQ,'logical',1,1,'storeVQ');
storeVtotQ = JJAsim_method_checkInput(storeVtotQ,'logical',1,1,'storeVtotQ');
storeEQ = JJAsim_method_checkInput(storeEQ,'logical',1,1,'storeEJtotQ');
repetitions = JJAsim_method_checkInput(repetitions,'double',1,1,'repetitions');
computePartitions = JJAsim_method_checkInput(computePartitions,'double',1,1,'computePartitions');
parallelQ = JJAsim_method_checkInput(parallelQ,'logical',1,1,'parallelQ');
cores = JJAsim_method_checkInput(cores,'double',1,1,'cores');

%get lists which entry into the input variables to use for each problem
switch inputMode
    case 'sweep'
        sweepW = [size(IExt,1),size(T,1),max([size(f,2),size(z,2)])];
        IExtProblemList = JJAsim_method_sweep(size(IExt,1),1,sweepW,repetitions);
        TProblemList = JJAsim_method_sweep(size(T,1),2,sweepW,repetitions);
        fProblemList = JJAsim_method_sweep(size(f,2),3,sweepW,repetitions);
        zProblemList = JJAsim_method_sweep(size(z,2),3,sweepW,repetitions);
    case 'list'
        IExtProblemList = JJAsim_method_list(size(IExt,1),W,repetitions);
        TProblemList = JJAsim_method_list(size(T,1),W,repetitions);
        fProblemList = JJAsim_method_list(size(f,2),W,repetitions);
        zProblemList = JJAsim_method_list(size(z,2),W,repetitions); 
end
th1ProblemList = JJAsim_method_list(size(th1,2),W,repetitions); 
if array.capacitanceQ
    th2ProblemList = JJAsim_method_list(size(th2,2),W,repetitions); 
end

%partition the problem lists in computePartitions partitions.
[IExtProblemSelection,IExtTableSelection,WIExt,partitionSize] = JJAsim_method_partition(IExtProblemList,computePartitions);
[TProblemSelection,TTableSelection,WT] = JJAsim_method_partition(TProblemList,computePartitions);
[fProblemSelection,fTableSelection,Wf] = JJAsim_method_partition(fProblemList,computePartitions);
[zProblemSelection,zTableSelection,Wz] = JJAsim_method_partition(zProblemList,computePartitions);
[th1ProblemSelection,th1TableSelection,Wth1] = JJAsim_method_partition(th1ProblemList,computePartitions);
if array.capacitanceQ
    [th2ProblemSelection,th2TableSelection,Wth2] = JJAsim_method_partition(th2ProblemList,computePartitions);
end


%get number of time points for output variables.
Ntth = sum(thTimePoints);    NtI = sum(ITimePoints);        NtE = sum(ETimePoints);             
NtV = sum(VTimePoints);      NtVtot = sum(VtotTimePoints);   
  

%initialize output variables
if storethQ;   	thOut = zeros(Nj,partitionSize,Ntth,computePartitions); end
if storeIQ;   	IOut = zeros(Nj,partitionSize,NtI,computePartitions); end
if storeVQ;  	VOut = zeros(Nj,partitionSize,NtV,computePartitions); end
if storeVtotQ; 	VtotOut = zeros(partitionSize,NtVtot,computePartitions); end
if storeEQ;     EOut = zeros(partitionSize,NtE,computePartitions); end

if cores == 0
    cores = feature('numcores');
end

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
        TP = JJAsim_method_tocell(T,TTableSelection(:,pInd),1,2);
        fP = JJAsim_method_tocell(f,fTableSelection(:,pInd),2,2);
        zP = JJAsim_method_tocell(z,zTableSelection(:,pInd),2,2);
        th1P = JJAsim_method_tocell(th1,th1TableSelection(:,pInd),2,2);
        IprobP = IExtProblemSelection(:,pInd);
        TprobP = TProblemSelection(:,pInd);
        fprobP = fProblemSelection(:,pInd);
        zprobP = zProblemSelection(:,pInd);
        th1probP = th1ProblemSelection(:,pInd);
        if array.capacitanceQ
            th2P = JJAsim_method_tocell(th2,th2TableSelection(:,pInd),2,2);
            th2probP = th2ProblemSelection(:,pInd);
        end   
        
        %prepare temporary variables    
        if storethQ;  	thOutTemp = zeros(size(thOut,1),size(thOut,2),size(thOut,3),cores); end
        if storeIQ;  	IOutTemp = zeros(size(IOut,1),size(IOut,2),size(IOut,3),cores); end
        if storeVQ;  	VOutTemp = zeros(size(VOut,1),size(VOut,2),size(VOut,3),cores); end
        if storeVtotQ;  VtotOutTemp = zeros(size(VtotOut,1),size(VtotOut,2),cores); end
        if storeEQ;     EOutTemp = zeros(size(EOut,1),size(EOut,2),cores); end
        
        %execute run. The algorithm depends on whether the array contains inductance and/or capacitance
        if ~array.inductanceQ && ~array.capacitanceQ
            parfor p = coreInd
                [thOutTemp(:,:,:,p),IOutTemp(:,:,:,p),VOutTemp(:,:,:,p),VtotOutTemp(:,:,p),EOutTemp(:,:,p)] = ...
                    JJAsim_2D_network_simulate_priv_noscr_nocap_CPU(...
                    Nn,Nj,Np,M,A,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,t,IExtP{p},IExtCompact,...
                    IprobP(:,p),WIExt(p),TP{p},TCompact,TprobP(:,p),WT,fP{p},fCompact,fprobP(:,p),Wf(p),zP{p},...
                    zCompact,zprobP(:,p),Wz(p),th1P{p},th1probP(:,p),Wth1(p),Ic,IcCompact,cpQ,cp,icp,...
                    Rn,RnCompact,storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,storeVQ,...
                    VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE);
            end
        end
        if array.inductanceQ && ~array.capacitanceQ
            parfor p = coreInd
                [thOutTemp(:,:,:,p),IOutTemp(:,:,:,p),VOutTemp(:,:,:,p),VtotOutTemp(:,:,p),EOutTemp(:,:,p)] = ...
                    JJAsim_2D_network_simulate_priv_scr_nocap_CPU(...
                    Nn,Nj,Np,M,A,LMode,betaL,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,t,IExtP{p},IExtCompact,...
                    IprobP(:,p),WIExt(p),TP{p},TCompact,TprobP(:,p),WT,fP{p},fCompact,fprobP(:,p),Wf(p),zP{p},...
                    zCompact,zprobP(:,p),Wz(p),th1P{p},th1probP(:,p),Wth1(p),Ic,IcCompact,cpQ,cp,icp,...
                    Rn,RnCompact,storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,...
                    storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,...
                    ETimePoints,NtE);
            end
        end
        if ~array.inductanceQ && array.capacitanceQ
            parfor p = coreInd
                [thOutTemp(:,:,:,p),IOutTemp(:,:,:,p),VOutTemp(:,:,:,p),VtotOutTemp(:,:,p),EOutTemp(:,:,p)] = ...
                    JJAsim_2D_network_simulate_priv_noscr_cap_CPU(...
                    Nn,Nj,Np,M,A,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,dt,IExtP{p},IExtCompact,IprobP(:,p),...
                    WIExt(p),TP{p},TCompact,TprobP(:,p),WT,fP{p},fCompact,fprobP(:,p),Wf(p),zP{p},zCompact,zprobP(:,p),...
                    Wz(p),th1P{p},th1probP(:,p),Wth1(p),th2P{p},th2probP(:,p),Wth2(p),Ic,IcCompact,cpQ,cp,icp,...
                    Rn,RnCompact,betaC,betaCCompact,Vscheme,storethQ,thTimePoints,Ntth,...
                    storeIQ,ITimePoints,NtI,...
                    storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE);
            end
        end
        if array.inductanceQ && array.capacitanceQ
            parfor p = coreInd
                [thOutTemp(:,:,:,p),IOutTemp(:,:,:,p),VOutTemp(:,:,:,p),VtotOutTemp(:,:,p),EOutTemp(:,:,p)] = ...
                    JJAsim_2D_network_simulate_priv_scr_cap_CPU(...
                    Nn,Nj,Np,M,A,LMode,betaL,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,dt,IExtP{p},IExtCompact,...
                    IprobP(:,p),WIExt(p),TP{p},TCompact,TprobP(:,p),WT,fP{p},fCompact,fprobP(:,p),Wf(p),zP{p},...
                    zCompact,zprobP(:,p),Wz(p),th1P{p},th1probP(:,p),Wth1(p),th2P{p},th2probP(:,p),Wth2(p),Ic,...
                    IcCompact,cpQ,cp,icp,Rn,RnCompact,betaC,betaCCompact,Vscheme,...
                    storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,...
                    storeVQ,VTimePoints,NtV,storeVtotQ,...
                    VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE);
            end
        end
        
        %store temporary variables to output variables
        if storethQ;  	thOut(:,:,:,pInd) = thOutTemp(:,:,:,coreInd); end
        if storeIQ;  	IOut(:,:,:,pInd) = IOutTemp(:,:,:,coreInd); end
        if storeVQ;  	VOut(:,:,:,pInd) = VOutTemp(:,:,:,coreInd); end
        if storeVtotQ;  VtotOut(:,:,pInd) = VtotOutTemp(:,:,coreInd); end
        if storeEQ;     EOut(:,:,pInd) = EOutTemp(:,:,coreInd); end    
    end
else    
    %execute simulation for all partitions on single core
    for p = 1:computePartitions
        
        %get data selection for this partition
        Iprob = IExtProblemSelection(:,p);    Isel = IExtTableSelection(:,p);
        Tprob = TProblemSelection(:,p);       Tsel = TTableSelection(:,p);
        fprob = fProblemSelection(:,p);       fsel = fTableSelection(:,p);
        zprob = zProblemSelection(:,p);       zsel = zTableSelection(:,p);
        th1prob = th1ProblemSelection(:,p);     th1sel = th1TableSelection(:,p);
        if array.capacitanceQ
            th2prob = th2ProblemSelection(:,p);     th2sel = th2TableSelection(:,p);
        end
        
        %execute simulation for this partition
        if ~array.inductanceQ && ~array.capacitanceQ
            [thOut(:,:,:,p),IOut(:,:,:,p),VOut(:,:,:,p),VtotOut(:,:,p),EOut(:,:,p)] = ...
                JJAsim_2D_network_simulate_priv_noscr_nocap_CPU(...
                Nn,Nj,Np,M,A,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,t,IExt(Isel,:),IExtCompact,...
                Iprob,WIExt(p),T(Tsel,:),TCompact,Tprob,WT,f(:,fsel),fCompact,fprob,Wf(p),z(:,zsel),...
                zCompact,zprob,Wz(p),th1(:,th1sel),th1prob,Wth1(p),Ic,IcCompact,cpQ,cp,icp,...
                Rn,RnCompact,storethQ,thTimePoints,Ntth,storeIQ,...
                ITimePoints,NtI,...
                storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE);
        end
        if array.inductanceQ && ~array.capacitanceQ
            [thOut(:,:,:,p),IOut(:,:,:,p),VOut(:,:,:,p),VtotOut(:,:,p),EOut(:,:,p)] = ...
                JJAsim_2D_network_simulate_priv_scr_nocap_CPU(...
                Nn,Nj,Np,M,A,LMode,betaL,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,t,IExt(Isel,:),...
                IExtCompact,Iprob,WIExt(p),T(Tsel,:),TCompact,Tprob,WT,f(:,fsel),fCompact,fprob,Wf(p),...
                z(:,zsel),zCompact,zprob,Wz(p),th1(:,th1sel),th1prob,Wth1(p),Ic,IcCompact,cpQ,cp,icp,...
                Rn,RnCompact,storethQ,thTimePoints,Ntth,...
                storeIQ,ITimePoints,NtI,...
                storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,...
                storeEQ,ETimePoints,NtE);
        end
        if ~array.inductanceQ && array.capacitanceQ
            [thOut(:,:,:,p),IOut(:,:,:,p),VOut(:,:,:,p),VtotOut(:,:,p),EOut(:,:,p)] = ...
                JJAsim_2D_network_simulate_priv_noscr_cap_CPU(...
                Nn,Nj,Np,M,A,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,dt,IExt(Isel,:),IExtCompact,...
                Iprob,WIExt(p),T(Tsel,:),TCompact,Tprob,WT,f(:,fsel),fCompact,fprob,Wf(p),z(:,zsel),...
                zCompact,zprob,Wz(p),th1(:,th1sel),th1prob,Wth1(p),th2(:,th2sel),th2prob,Wth2(p),Ic,...
                IcCompact,cpQ,cp,icp,Rn,RnCompact,betaC,betaCCompact,Vscheme,...
                storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,...
                storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,...
                NtVtot,storeEQ,ETimePoints,NtE);
        end
        if array.inductanceQ && array.capacitanceQ
            [thOut(:,:,:,p),IOut(:,:,:,p),VOut(:,:,:,p),VtotOut(:,:,p),EOut(:,:,p)] = ...
                JJAsim_2D_network_simulate_priv_scr_cap_CPU(...
                Nn,Nj,Np,M,A,LMode,betaL,partitionSize,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,dt,IExt(Isel,:),...
                IExtCompact,Iprob,WIExt(p),T(Tsel,:),TCompact,Tprob,WT,f(:,fsel),fCompact,fprob,Wf(p),...
                z(:,zsel),zCompact,zprob,Wz(p),th1(:,th1sel),th1prob,Wth1(p),th2(:,th2sel),th2prob,Wth2(p),...
                Ic,IcCompact,cpQ,cp,icp,Rn,RnCompact,betaC,betaCCompact,Vscheme,...
                storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,...
                storeVQ,VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE);
        end
    end
end

%store output variables in out 
if storethQ
    thOut = permute(thOut,[1,2,4,3]);
    thOut = reshape(thOut,Nj,partitionSize*computePartitions,Ntth);
    out.th = thOut(:,1:Wtot,:); 
end
if storeIQ
    IOut = permute(IOut,[1,2,4,3]);
    IOut = reshape(IOut,Nj,partitionSize*computePartitions,NtI);
    out.I = IOut(:,1:Wtot,:); 
end
if storeVQ
    VOut = permute(VOut,[1,2,4,3]);
    VOut = reshape(VOut,Nj,partitionSize*computePartitions,NtV);
    out.V = VOut(:,1:Wtot,:); 
end
if storeVtotQ
    VtotOut = permute(VtotOut,[1,3,2]);
    VtotOut = reshape(VtotOut,partitionSize*computePartitions,NtVtot);
    out.Vtot = VtotOut(1:Wtot,:); 
end
if storeEQ
    EOut = permute(EOut,[1,3,2]);
    EOut = reshape(EOut,partitionSize*computePartitions,NtE);
    out.E = EOut(1:Wtot,:); 
end
end
