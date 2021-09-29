function [thOut,nOut,IOut,EOut,solutionQ,nrOfSteps] = ...
    JJAsim_2D_network_stationairyState_priv_noscr_CPU(Nn,Nj,Np,M,A,...
    W,IExtBaseJ,pathArea,areaCompact,IExt,IExtprob,WIExt,f,fCompact,fprob,Wf,...
    z,zCompact,zprob,Wz,ic,icprob,Wic,Ic,IcCompact,cpQ,cp,dcp,icp,...
    storeThQ,storenQ,storeIQ,storeEQ,newton_steps,newton_tol) %#ok<*INUSL>

%input:
% Nn               1 by 1            	number of nodes
% Nj                1 by 1            	number of junctions
% Np                1 by 1             	number of paths
% M                 Nn by Nj        	Kirchhoffs current rules
% A                 Np by Nj           	Kirchhoffs voltage rules
% W                 1 by 1            	number of parallel problems
% IExtBaseJ         Nj by 1           	Junction external current base
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
% cpQ               1 by 1              true if custom current phase relation
% cp                function_handle     custom current phase relation
% dcp               function_handle     th-derivative of custom cp relation
% icp               function_handle     th-integral of custom cp relation
% storeThQ          1 by 1          	if true, store th data
% storenQ           1 by 1           	if true, store n data
% storeIQ           1 by 1           	if true, store I data
% storeEQ           1 by 1            	if true, store E data
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

%find stationairy states
%using newton method aiming to solve the system f(th) = 0 with:
%f = [
%  M*(cp(Ic,th) - Ip)
%  A*th - 2*pi*(z - f)
%    ];

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
        Q = (cp(Ic,th) - Ip)./P;
    else
        Q = (Ic.*sin(th) - Ip)./P;
    end
    
    J = zeros(Np,size(th,2));
    DP1=[];
    for w = 1:size(th,2)
        size(A);
        sizeP=size(1./P(:,w));
        sizeAT=size(AT);
        AT;
        [ii,jj,ss]=find(AT);
        sizess=size(ss);
        P1=1./P(:,w);
        for kk = 1:sizess(1)
            d = ii(kk);
            ss(kk) = ss(kk)*P1(d);
        end
        
%         1./P(:,w)
%         1./P(:,w).*AT
%         ss = (1./P(:,w))*ss;
        DP=sparse(ii,jj,ss);
        DP=A*DP;
        
%         DP = 1./P(:,w).*AT;
%         size(DP)
%         
%         DP = times(oneover,AT);
%         size(DP);
%         DP = A*DP;
%         for 
%         size(P1)
%         size(AT)
%         %producat = P1.*AT
%         sizeAT = size(AT);
%         for k=1:sizeAT(2)
%         P1=1./P(:,w);
%         P2=P1(k);
%         DP1(:,k) = P2.*AT(:,k);
%         end
%         DP=A*(1./P(:,w).*AT);
        %sizeDP=size(DP)
        %sizeAQ=size(A*Q(:,w))
        J(:,w) = DP\(A*Q(:,w));
    end
    th = th + ((AT*J)./P - Q);
    if cpQ
        notConverged = mean(abs(M*(cp(Ic,th) - Ip)),1) > newton_tol;
    else
        notConverged = mean(abs(M*(Ic.*sin(th) - Ip)),1) > newton_tol;
    end
    %notConverged = notConverged | (mean(abs(mod(A*th - df,2*pi))) > newton_tol);
    selnp(:,problemSelection) = ~notConverged;
    thOut(:,selnp) = th(:,~notConverged);
    th = th(:,notConverged);
    Ip = Ip(:,notConverged);
    selp(:,problemSelection) = notConverged;
    problemSelection = problemSelection & selp;
    newtonStep = newtonStep + 1;
    nrOfSteps(problemSelection) = nrOfSteps(problemSelection) + 1;
end
thOut(:,problemSelection) = th;

%get list for which problems a solution is found. 
solutionQ = ~problemSelection;

%assign store output variables
if storenQ
    nOut = z(:,zprob)-A*round(thOut/2/pi);
else
    nOut = [];
end
if storeIQ   
    if cpQ
        IOut = cp(Ic,thOut);       
    else
        IOut = Ic.*sin(thOut);                 
    end
else
    IOut = []; 
end
if storeEQ
    if cpQ
        EOut = sum(icp(Ic,thOut),1)'; 
    else
        EOut = sum(Ic.*(1-cos(thOut)),1)'; 
    end
else
    EOut = []; 
end
if ~storeThQ   
    thOut = []; 
end
end