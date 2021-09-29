function [thOut,nOut,IOut,EOut,solutionQ,nrOfSteps] = ...
    JJAsim_2D_network_stationairyState_priv_scr_CPU(Nn,Nj,Np,M,A,betaL,LMode,...
    W,IExtBaseJ,pathArea,areaCompact,IExt,IExtprob,WIExt,f,fCompact,fprob,Wf,...
    z,zCompact,zprob,Wz,ic,icprob,Wic,Ic,IcCompact,cpQ,cp,dcp,icp,storethQ,...
    storenQ,storeIQ,storeEQ,newton_steps,newton_tol) %#ok<*INUSL>
%
%using newton method aiming to solve the system f(th) = 0 with:
%f = [
%  M*cp(Ic,th) - IExt
%  A*th - 2*pi*(z - f) - A*betaL*Icp(Ic,th))
%    ];
%
%input:
% Nn               1 by 1            	number of nodes
% Nj                1 by 1            	number of junctions
% Np                1 by 1             	number of paths
% M                 Nn by Nj        	Kirchhoffs current rules
% A                 Np by Nj           	Kirchhoffs voltage rules
% betaL             1 by 1           	Inductance, if LMode = 0, uniform self inductance
%                   or Nj by Nj                     if LMode = 1, sparse/dense inductance matrix
% LMode             1 by 1              0 or 1, determines how inductance is specified.
% W                 1 by 1            	number of parallel problems
% IExtBaseJ         Nj by 1           	junction external current base
% pathArea          (Np or 1) by 1    	area enclosed by paths
% areaCompact       1 by 1           	shows which dimensions of pathArea are compact 
%                                      	(i.e. singleton because it is otherwise the same value repeated)
% IExt              WIExt by 1          external current table
% Iprob             W by 1           	List for each problem which row from IExt is taken
% WIExt             1 by 1            	number of rows in the IExt table    
% f                 (Np or 1) by Wf    	frustration factor table
% fCompact          1 by 2            	which dimensions of f are compact
% fprob             W by 1             	List for each problem which column from f is taken
% Wf                1 by 1             	number of columns in the f table   
% z                 (Np or 1) by Wz    	phase zone table
% zCompact          1 by 2             	which dimensions of z are compact
% zProb             W by 1            	List for each problem which column from z is taken
% Wz                1 by 1            	number of columns in the z table   
% ic                Nj by Wic       	initial phase difference table
% icprob            W by 1            	List for each problem which column from ic is taken
% Wic               1 by 1            	number of columns in the ic table            
% Ic                (Nj or 1) by 1     	critical current values
% IcCompact         1 by 1             	which dimensions of Ic are compact  
% storethQ          1 by 1          	if true, store th data
% storenQ           1 by 1           	if true, store n data
% storeIQ           1 by 1           	if true, store I data
% storeEQ           1 by 1            	if true, store EJtot data
% newton_steps      1 by 1              number of steps of newton algorithm
% newton_tol        1 by 1              acceptance tolerance of steps of newton algorithm

%prepare quantities
AT = A';

%get initial condition
df = 2*pi*(z(:,zprob)-pathArea.*f(:,fprob));
dfCompact = fCompact(1) && zCompact(1)  && areaCompact;
if dfCompact
    df = repmat(df,Np,1);
end
th = ic(:,icprob);

%get external current
Ip =  shiftdim(IExt(IExtprob),-1).*IExtBaseJ;

%initialize stationairy state finding
selp = true(1,W);
nrOfSteps = ones(1,W);
thOut = zeros(size(th));
newtonStep = 1;
problemSelection = true(1,W);
J0 = df - A*(betaL*Ip);
dfp = df;

%resulting in the iterative scheme: 
while sum(problemSelection) > 0 && newtonStep <= newton_steps
    selnp = false(1,W);

    if cpQ
        P = dcp(Ic,th);
    else
        P = Ic.*cos(th);
    end
    P(abs(P) < 1E-10) = P(abs(P) < 1E-10) + 1E-2;
    if cpQ
        Q = (cp(Ic,th) - Ip);
    else
        Q = (Ic.*sin(th) - Ip);
    end
    
    J = J0 + A*(Q./P - th);
    for w = 1:size(th,2)
        if LMode == 0
            DP = A*((1./P(:,w) + betaL).*AT);
        else
            DP = A*((spdiags(1./P(:,w),0,Nj,Nj) + betaL)*AT);
        end
        J(:,w) = DP\J(:,w);
    end
    th = th + ((AT*J- Q)./P);
    if cpQ
        I = cp(Ic,th);
    else
        I = Ic.*sin(th);
    end
    notConverged = mean(abs(M*(I - Ip)),1) > newton_tol;
    notConverged = notConverged | (mean(abs(A*(th + betaL*I) - dfp)) > newton_tol);
    selnp(:,problemSelection) = ~notConverged;
    thOut(:,selnp) = th(:,~notConverged);
    th = th(:,notConverged);
    Ip = Ip(:,notConverged);
    J0 = J0(:,notConverged);
    dfp = dfp(:,notConverged);
    selp(:,problemSelection) = notConverged;
    problemSelection = problemSelection & selp;
    newtonStep = newtonStep + 1;
    nrOfSteps(problemSelection) = nrOfSteps(problemSelection) + 1;
end
thOut(:,problemSelection) = th;

%get list for which problems a solution is found. 
solutionQ = ~problemSelection;

%compute junction current if necessary
if storeIQ || storeEQ
    if cpQ
        I = cp(Ic,thOut);
    else
        I = Ic.*sin(thOut);
    end
end

%store vortex configuration
if storenQ
    nOut = z(:,zprob)-A*round(thOut/2/pi);
else
    nOut = []; 
end

%store junction current
if storeIQ   
    IOut = I;      
else
    IOut = []; 
end

%store energy
if storeEQ
    if cpQ
        EOut = sum(icp(Ic,thOut),1)'; 
    else
        EOut = sum(Ic.*(1-cos(thOut)),1)'; 
    end
    if LMode == 0
        EOut =EOut + 0.5*betaL*sum(I.^2,1)'; 
    else
        EOut =EOut + 0.5*sum(I.*(betaL*I),1)'; 
    end
else
    EOut = []; 
end

%store th
if ~storethQ
    thOut = []; 
end
end

