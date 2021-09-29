function out = JJAsim_2D_network_method_getIExtMax(array,n,f,varargin)
%out = JJAsim_2D_network_method_getIExtMax(array,n,f,varargin)
%
%DESCRIPTION
% - Finds the maximal external current that can be sourced through an array
%   without resistance.
% - Starts at IExt = 0, and increses it every step with stepSize, untill
%   it fails to find a solution at IExt. Then, stepSize is halved and 
%   the process continues, for a total number of steps of nrOfSteps.
%   This process will result in a lower and upper bound.
% - The output contains IExtLowerBound and IExtUpperBound. The output
%   flag determines the status of the solution (see below).
% - The output th is the stationairy state solution at IExtLowerBound.
%   Note, this is given in the phaze zone z = 0.
% - One can manually give an initial guess with th1. Note this is also
%   in the phase zone z = 0.
% - Several problems can be computed in parallel in one function call.
% - All physical quantities are normalized, see documentation.
%
% flag  0:  Solution found between IExtLowerBound and IExtUpperBound
%       1:  Solution found but no upperbound. Consider increasing nrOfSteps 
%           or stepSize
%       2:  No solution found, even for IExt = 0. Consider manual th1.
%
%FIXED INPUT
% array                 struct      Geometric information of the network. 
%                                   Create with JJAsim_2D_network_create.       
% n                     Np by W*  	Phase zone
% f                     Np* by W* 	Frustration factor
%
%VARIABLE INPUT
% initialGuessAlgorithm string      'arctan' or 'london' or 'manual'
% th1                   Nj by W    	Manual initial guess for IExt = 0. In phase 
%                                   zone z = 0.
% nrOfSteps             1 by 1    	number of steps.
% stepSize              1 by 1    	Initial step size. Gets refined when upper 
%                                   bound is found.
% newton_steps          1 by 1      number of steps for newtons algorithm.
% newton_tol            1 by 1      acceptance tolerance for newtons algorithm.
% computeAlgorithm      string      'static' or 'dynamic'. How is determined if at
%                                   an IExt value a solution exists. For 'static'
%                                   this is done by using the stationairyState 
%                                   function, for 'dynamic' this is done by a time
%                                   evolution, and checking if the second half is
%                                   unchanged. In both cases it is also required that 
%                                   the vortex configuration is equal to n.
%  (if static)
%     checkStabilityQ   1 by 1      If true, a state must also be dynamically stable
%                                   to count as an existing solution (slow).
%  (if dynamic)
%     t                 Nt by 1     Time points that are used for time evolution.
%     intervalSize      1 by 1      Size (between 0 and 1) of time interval used to 
%                                   look if th is unchanged
%     thTol                         maximum deviation th can have in the interval
% computePartitions     1 by 1     	Splits the problem set into computePartitions 
%                                   which are computed consecutively. Use to reduce 
%                                   memory usage or for parallel computation. 
% parallelQ             1 by 1    	If true, the partitions are computed in 
%                                   parallel on multiple cores.  
% cores                 1 by 1    	Number of cores used by parallel computing. If 0, 
%                                   it is set to the amount of available cores.            
%
%OUTPUT
% out.flag              1 by W    	true if valid solution is found
% out.IExtLowerBound    1 by W      The maximum IExt is atleast IExtLowerBound
% out.IExtUpperBound    1 by W      The maximum IExt < IExtUpperBound
% out.th                Nj by W    	phase configuration of solution at IExt = 
%                                   IExtLowerBound in phase zone z = 0.

%get variable input
inputParameters = {
    'initialGuessAlgorithm'  'arctan'
    'th1'                    []    
    'nrOfSteps'              30
    'stepSize'               0.4
    'newton_steps'         	 20
    'newton_tol'           	 1E-10
    'computeAlgorithm'       'static'
    'checkStabilityQ'        true
    't'                      []
    'intervalSize'           0.5   
    'thTol'                  0.01
    'computePartitions'      1
    'parallelQ'              false
    'cores'                  0
    'Istart'                  0
    };

options = JJAsim_method_parseOptions(inputParameters,varargin,...
    'JJAsim_2D_network_stationairyState');
initialGuessAlgorithm = options.initialGuessAlgorithm;
th1 = options.th1;
nrOfSteps = options.nrOfSteps;
stepSize = options.stepSize;
newton_steps = options.newton_steps;
newton_tol = options.newton_tol;
computeAlgorithm = options.computeAlgorithm;
checkStabilityQ = options.checkStabilityQ;
t = options.t;
intervalSize = options.intervalSize;
thTol = options.thTol;
computePartitions = options.computePartitions;
parallelQ = options.parallelQ;
cores = options.cores;
Istart = options.Istart;

%get sizes
Nj = array.Nj;
Np = array.Np;

%check if the array struct corresponds to that of a 2D network with free boundaries.
if ~(strcmp(array.type,'network') && array.ndims == 2)
   error('array not 2D network'); 
end

%get number of problems
W = max([size(f,2),size(n,2)]);

%check input size and type
f = JJAsim_method_checkInput(f,'double',[Np,W],[1,1],'f');
n = JJAsim_method_checkInput(n,'double',[Np,W],[0,1],'n');
if size(f,1) == 1
    f = repmat(f,Np,1);
end
if size(f,2) == 1
    f = repmat(f,1,W);
end
if size(n,2) == 1
    n = repmat(n,1,W);
end

%check other input
nrOfSteps = JJAsim_method_checkInput(nrOfSteps,'double',1,1,'nrOfSteps');
stepSize = JJAsim_method_checkInput(stepSize,'double',1,1,'stepSize');
computePartitions = JJAsim_method_checkInput(computePartitions,'double',1,1,...
    'computePartitions');
parallelQ = JJAsim_method_checkInput(parallelQ,'logical',1,1,'parallelQ');
cores = JJAsim_method_checkInput(cores,'double',1,1,'cores');
if cores == 0
    cores = feature('numcores');
end
switch computeAlgorithm
    case 'dynamic'
        Nt = length(t);
        t = JJAsim_method_checkInput(t,'double',Nt,0,'t');
        tP = round(Nt*(1-intervalSize));
end


%set parallel pool
if parallelQ
    %initialize the parallel pool
    JJAsim_method_poolinit(cores);
else
    %remove any present parallel pool to make up free space.
    JJAsim_method_poolinit(1);
end

%get initial guess
switch initialGuessAlgorithm
    case 'manual'
        th1 = JJAsim_method_checkInput(th1,'double',[Nj,W],[0,0],'th1');

    case 'arctan'
        L = max(sum(n~=0,1));
        x0 = zeros(L,W); y0 = zeros(L,W); n0 = zeros(L,W);
        for i = 1:W
            ind = n(:,i) ~= 0;
            x0(1:sum(ind),i) =  array.pathCentroid(ind,1);
            y0(1:sum(ind),i) =  array.pathCentroid(ind,2);
            n0(1:sum(ind),i) =  n(ind,i);
        end
        th1 = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f);
        
    case 'london'
        th1 = JJAsim_2D_network_stationairyState_approx_london(array,n,f);
        th1 = JJAsim_2D_network_method_changePhaseZone(array,th1,n,0);
        
    otherwise
        error('unrecognized initial guess algorithm');
end

%array.A*th1
%get initial exact states
switch computeAlgorithm
    case 'static'
        out = JJAsim_2D_network_stationairyState(array,'list',Istart,f,0,th1,'computePartitions',...
            computePartitions,'parallelQ',parallelQ,'cores',cores,'newton_steps',newton_steps,...
            'newton_tol',newton_tol);
        nOut = JJAsim_2D_network_method_getn(array,out.th,0);
        if checkStabilityQ
            stableQ = JJAsim_2D_network_stationairyState_stability(array,out.th);
            solQ = out.solutionQ' & (mean(abs(nOut-n),1) < 1E-10) & stableQ;
        else
            solQ = out.solutionQ' & (mean(abs(nOut-n),1) < 1E-10);
        end
        th1 = out.th;
        
    case 'dynamic'
        out = JJAsim_2D_network_simulate(array,t,'list',0,f,0,th1,'computePartitions',...
            computePartitions,'parallelQ',parallelQ,'cores',cores);
        nOut = JJAsim_2D_network_method_getn(array,out.th(:,:,end),0);
        thi = out.th(:,:,tP:end);
        crit1 = sum(sum(thi -min(thi,3) > thTol/2,3),1) ~= 0;
        crit2 = sum(sum(max(thi,3) -thi > thTol/2,3),1) ~= 0;
        solQ = crit1 & crit2 & (mean(abs(nOut-n),1) < 1E-10);
        th1 = thi(:,:,end);
end

%handle initial guesses that did not result in a solution even at IExt = 0. 
% by setting flag = 2 and removing them from consideration.
flag = nan(1,W);
flag(solQ == false) = 2;
selection = solQ;
f = f(:,selection);
n = n(:,selection);
th1 = th1(:,selection);
Wp = sum(selection);

clear out;

if Wp ~= 0

%start IExt finding
step = ones(1,Wp)*stepSize;
upperboundFoundQ = false(1,Wp);
IExt = Istart + step;
outSteps = zeros(nrOfSteps,Wp);
outIExt = zeros(nrOfSteps,Wp);
for i = 1:nrOfSteps
    switch computeAlgorithm
        case 'static'
            outIExt(i,:) = IExt;
            out = JJAsim_2D_network_stationairyState(array,'list',IExt',f,0,th1,...
                'computePartitions',computePartitions,'parallelQ',parallelQ,'cores',...
                cores,'newton_steps',newton_steps,'newton_tol',newton_tol);
            outSteps(i,:) = out.nrOfSteps';
            nOut = JJAsim_2D_network_method_getn(array,out.th,0);
            if checkStabilityQ
                stableQ = JJAsim_2D_network_stationairyState_stability(array,out.th);
                solQ = out.solutionQ' & (mean(abs(nOut-n),1) < 1E-10) & stableQ;
            else
                solQ = out.solutionQ' & (mean(abs(nOut-n),1) < 1E-10);
            end
            th1(:,solQ) = out.th(:,solQ);
            
        case 'dynamic'
            out = JJAsim_2D_network_simulate(array,t,'list',IExt',f,0,th1,'computePartitions',...
                computePartitions,'parallelQ',parallelQ,'cores',cores);
            nOut = JJAsim_2D_network_method_getn(array,out.th(:,:,end),0);
            thi = out.th(:,:,tP:end);
            crit1 = sum(sum(thi -min(thi,3) > thTol/2,3),1) ~= 0;
            crit2 = sum(sum(max(thi,3) -thi > thTol/2,3),1) ~= 0;
            solQ = crit1 & crit2 & (mean(abs(nOut-n),1) < 1E-10);
            th1 = thi(:,:,end);
    end
    IExt(solQ) = IExt(solQ) + step(solQ);
    IExt(~solQ) = IExt(~solQ) - step(~solQ)/2;
    step(~solQ) = step(~solQ)/2;
    upperboundFoundQ = upperboundFoundQ | (~solQ);
end
IExtLowerBound = IExt - step;
IExtUpperBound = IExt;
flagp = zeros(1,Wp);
flagp(~upperboundFoundQ) = 1;
IExtUpperBound(~upperboundFoundQ) = inf;

%assign output data
clear out;

out.IExtLowerBound = zeros(1,W);
out.IExtLowerBound(selection) = IExtLowerBound;
out.IExtUpperBound = zeros(1,W);
out.IExtUpperBound(selection) = IExtUpperBound;
out.steps = outSteps;
out.IExtEvolution = outIExt;
flag(selection) = flagp;
end

out.flag = flag;
out.th = nan(Nj,W);
out.th(:,selection) = th1;

end