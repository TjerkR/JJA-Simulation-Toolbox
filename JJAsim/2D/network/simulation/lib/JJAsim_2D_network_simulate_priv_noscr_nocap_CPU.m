function [thOut,IOut,VOut,VtotOut,EOut] = JJAsim_2D_network_simulate_priv_noscr_nocap_CPU(...
    Nn,Nj,Np,M,A,W,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,t,IExt,IExtCompact,...
    IExtprob,WIExt,T,TCompact,Tprob,WT,f,fCompact,fprob,Wf,z,zCompact,zprob,Wz,ic,icprob,Wic,Ic,...
    IcCompact,cpQ,cp,icp,Rn,RnCompact,storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,storeVQ,...
    VTimePoints,NtV,storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE) %#ok<*INUSL>

%input:
% Nn               1 by 1            	number of nodes
% Nj                1 by 1            	number of junctions
% Np                1 by 1             	number of paths
% M                 Nn by Nj        	Kirchhoffs current rules
% A                 Np by Nj           	Kirchhoffs voltage rules
% W                 1 by 1            	number of parallel problems
% IExtBaseJ         Nj by 1             Junction external current base
% IExtBaseTotal     1 by 1              Total unit injected current with IExtBase
% pathArea          (Np or 1) by 1    	area enclosed by paths
% areaCompact       1 by 1           	shows which dimensions of pathArea are compact 
%                                      	(i.e. singleton because it is otherwise the same value repeated)
% Nt                1 by 1            	number of time points
% t                 Nt by 1           	time values of all time points
% IExt              WIExt by (Nt or 1)	external current table
% IExtCompact       1 by 2            	which dimensions of IExt are compact
% IExtprob          W by 1           	List for each problem which row from IExt is taken
% WIExt             1 by 1            	number of rows in the IExt table    
% T                 WT by (Nt or 1)  	Temperature table
% TCompact          1 by 2             	which dimensions of T are compact
% Tprob             W by 1             	List for each problem which row from T is taken
% WT                1 by 1             	number of rows in the T table   
% f                 (Np or 1) by Wf    	frustration factor table
% fCompact          1 by 2            	which dimensions of f are compact
% fprob             W by 1             	List for each problem which column from f is taken
% Wf                1 by 1             	number of columns in the f table   
% z                 (Np or 1) by Wz    	phase zone table
% zCompact          1 by 2             	which dimensions of z are compact
% zProb             W by 1            	List for each problem which column from z is taken
% Wz                1 by 1            	number of columns in the z table   
% ic                Nj by Wic       	initial phase differences table 
% icprob            W by 1            	List for each problem which column from ic table is taken
% Wic               1 by 1            	number of columns in the ic table       
% Ic                (Nj or 1) by 1     	critical current values
% IcCompact         1 by 1             	which dimensions of Ic are compact
% cpQ               1 by 1              true if custom current phase relation
% cp                function_handle     custom current phase relation
% icp               function_handle     integral of custom current phase relation
% Rn                (Nj or 1) by 1   	normal state resistance values
% RnCompact         1 by 1            	which dimensions of Rn are compact
% storethQ          1 by 1          	if true, store th data
% thTimePoints      Nt by 1          	which timepoints are stored in the th output
% Ntth              1 by 1            	number of timepoints in th output
% storenQ           1 by 1           	if true, store n data
% nTimePoints       Nt by 1           	which timepoints are stored in the n output
% Ntn               1 by 1            	number of timepoints in n output
% storeIQ           1 by 1           	if true, store I data
% ITimePoints       Nt by 1             which timepoints are stored in the I output
% NtI               1 by 1           	number of timepoints in U output
% storeVQ           1 by 1            	if true, store V data
% VTimePoints       Nt by 1           	which timepoints are stored in the V output
% NtV               1 by 1          	number of timepoints in V output
% storeVtotQ     	1 by 1           	if true, store Vtot data
% VtotTimePoints  	Nt by 1          	which timepoints are stored in the Vtot output
% NtVtot            1 by 1          	number of timepoints in Vtot output
% storeEQ           1 by 1            	if true, store E data
% ETimePoints       Nt by 1           	which timepoints are stored in the E output
% NtE               1 by 1             	number of timepoints in E output

%initialize output variables
if storethQ;    thOut = zeros(Nj,W,Ntth);    else; thOut = []; end
if storeIQ;     IOut = zeros(Nj,W,NtI);      else; IOut = []; end
if storeVQ;     VOut = zeros(Nj,W,NtV);      else; VOut = []; end
if storeVtotQ;  VtotOut = zeros(W,NtVtot);   else; VtotOut = []; end
if storeEQ;     EOut = zeros(W,NtE);         else; EOut = []; end

thTimePoint = 1;       ETimePoint = 1;
ITimePoint = 1;        VTimePoint = 1;     
VtotTimePoint = 1;

%compute matrices and geometric quantities
AT = A';
D = A*AT;
DR = A*(Rn.*AT);
df = 2*pi*(z(:,zprob)-pathArea.*f(:,fprob));
dfCompact = fCompact(1) && zCompact(1) && areaCompact(1);
if dfCompact
    g = AT*(D\repmat(df,Np,1));
else
    g = AT*(D\df);
end

%handle initial condition
thnp = ic(:,icprob);

%shift dimension of T and IExt to the right by 1 to line up with the rest
T = shiftdim(T,-1);
IExt = shiftdim(IExt,-1);

%do computations 
for i = 1:Nt-1
    dt = t(i+1)-t(i);
    thn = thnp;

    %correct th every 100 steps
    if mod(i-1,100) == 0
        thn = thn - AT*(D\(A*(thn-g)));
    end
    
    %check if I and Temp are time dependent
    if ~TCompact(2); iTemp = i; else; iTemp = 1; end
    if ~IExtCompact(2); iI = i; else; iI = 1; end

    %store th and E for this timestep
    if storethQ && thTimePoints(i)
        thOut(:,:,thTimePoint) = thn;
        thTimePoint = thTimePoint + 1;
    end
    if storeEQ && ETimePoints(i)
        if cpQ
            EOut(:,ETimePoint) = shiftdim(sum(icp(Ic,thn),1),1);
        else
            EOut(:,ETimePoint) = shiftdim(sum(Ic.*(1-cos(thn)),1),1);
        end
        ETimePoint = ETimePoint + 1;
    end
    
    %core computations of this timestep
    IEn = IExtBaseJ.*IExt(1,IExtprob,iI); 
    if cpQ
        yn = cp(Ic,thn);
    else
        yn = Ic.*sin(thn);
    end
    yn = yn + randn(Nj,W).*sqrt(2*T(1,Tprob,iTemp)/dt./Rn)- IEn;
    thnp = AT*(DR\(A*(Rn.*yn)));
    thnp = thn - dt*Rn.*(yn-thnp);
    
    %store I, V and Vtot for this timestep
    if (storeIQ && ITimePoints(i)) || (storeVQ && VTimePoints(i)) || (storeVtotQ && VtotTimePoints(i))
        Vn = (thnp-thn)/dt;
    end
    if storeIQ && ITimePoints(i)
        IOut(:,:,ITimePoint) = yn + IEn + Vn./Rn;
        ITimePoint = ITimePoint + 1;
    end
    if storeVQ && VTimePoints(i)
        VOut(:,:,VTimePoint) = Vn;
        VTimePoint = VTimePoint + 1;
    end
    if storeVtotQ && VtotTimePoints(i)
        VtotOut(:,VtotTimePoint) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
        VtotTimePoint = VtotTimePoint + 1;
    end
end

%store th and E for last timestep
if storethQ && thTimePoints(Nt)
    thOut(:,:,end) = thnp;
end
if storeEQ && ETimePoints(Nt)
    if cpQ
        EOut(:,end) = shiftdim(sum(icp(Ic,thnp),1),1);
    else
        EOut(:,end) = shiftdim(sum(Ic.*(1-cos(thnp)),1),1);
    end
end
end

