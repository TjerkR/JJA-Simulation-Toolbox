function [thOut,IOut,VOut,VtotOut,EOut] = JJAsim_2D_network_simulate_priv_scr_cap_CPU(...
    Nn,Nj,Np,M,A,LMode,betaL,W,IExtBaseJ,IExtBaseTotal,pathArea,areaCompact,Nt,dt,IExt,...
    IExtCompact,IExtprob,WIExt,T,TCompact,Tprob,WT,f,fCompact,fprob,Wf,z,zCompact,zprob,Wz,...
    ic1,ic1prob,Wic1,ic2,ic2prob,Wic2,Ic,IcCompact,cpQ,cp,icp,Rn,RnCompact,betaC,betaCCompact,...
    scheme,storethQ,thTimePoints,Ntth,storeIQ,ITimePoints,NtI,storeVQ,VTimePoints,NtV,...
    storeVtotQ,VtotTimePoints,NtVtot,storeEQ,ETimePoints,NtE) %#ok<*INUSL>

%input:
% Nn               1 by 1            	number of nodes
% Nj                1 by 1            	number of junctions
% Np                1 by 1             	number of paths
% M                 Nn by Nj        	Kirchhoffs current rules
% A                 Np by Nj           	Kirchhoffs voltage rules
% LMode             1 by 1              inductance mode. determines size of betaL
% betaL (LMode = 0) 1 by 1              constant self inductance (all junctions have same betaL)
%       (LMode = 1) Nj by Nj            junction-inductance matrix. sparse or dense.
% W                 1 by 1            	number of parallel problems
% IExtBaseJ         Nj by 1             Junction external current base
% IExtBaseTotal     1 by 1              Total unit injected current with IExtBase
% pathArea          (Np or 1) by 1    	area enclosed by paths
% areaCompact       1 by 1           	shows which dimensions of pathArea are compact 
%                                      	(i.e. singleton, implying it is the same value repeated)
% Nt                1 by 1            	number of time points
% dt                1 by 1           	timestep
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
% ic1               Nj by Wic1       	timestep 1 phase differences table
% ic1prob           W by 1            	List for each problem which column from ic1 table is taken
% Wic1              1 by 1            	number of columns in the ic1 table    
% ic2               Nj by Wic2      	timestep 2 phase differences table
% ic2prob           W by 1            	List for each problem which column from ic2 table is taken
% Wic2              1 by 1            	number of columns in the ic2 table
% Ic                (Nj or 1) by 1     	critical current values
% IcCompact         1 by 1             	which dimensions of Ic are compact
% cpQ               1 by 1              true if custom current phase relation
% cp                function_handle     custom current phase relation
% icp               function_handle     integral of custom current phase relation
% Rn                (Nj or 1) by 1   	normal state resistance values
% RnCompact         1 by 1            	which dimensions of Rn are compact
% betaC             (Nj or 1) by 1   	mccumber beta parameter for capacitance
% betaCCompact      1 by 1            	which dimensions of betaC are compact
% scheme            1 by 1              dth/dt scheme; 0 for forward, 1 for central and 2 for backward
% storethQ          1 by 1          	if true, store th data
% thTimePoints      Nt by 1          	which timepoints are stored in the th output
% Ntth              1 by 1            	number of timepoints in th output
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

thTimePoint = 1;    ETimePoint = 1;     ITimePoint = 1;     
VTimePoint = 1;     VtotTimePoint = 1; 

%compute matrices and geometric quantities
AT = A';
D = A*AT;

%get inductance related matrices
if LMode == 1
    DB = A*betaL*AT;
end

%get resistance and capacitance matrices
switch scheme
    case 0
        Cp = betaC/dt^2 + 1/dt./Rn;
        C0 = -2*betaC/dt^2 - 1/dt./Rn;
        Cm = betaC/dt^2;
    case 1
        Cp = betaC/dt^2 + 1/dt/2./Rn;
        C0 = -2*betaC/dt^2;
        Cm = betaC/dt^2 - 1/dt/2./Rn;
    case 2
        Cp = betaC/dt^2;
        C0 = -2*betaC/dt^2 + 1/dt./Rn;
        Cm = betaC/dt^2 - 1/dt./Rn;
end

%get gauge
df = 2*pi*(z(:,zprob)-pathArea.*f(:,fprob));
dfCompact = fCompact(1) && zCompact(1) && areaCompact(1);
if dfCompact
    g = AT*(D\repmat(df,Np,1));
else
    g = AT*(D\df);
end

%get total injected current
ILBase = betaL*IExtBaseJ;

%handle initial condition
thm = ic1(:,ic1prob);
th0 = ic2(:,ic2prob);

%shift dimension of T and IExt to the right by 1 to line up with the rest
T = shiftdim(T,-1);
IExt = shiftdim(IExt,-1);

%store instantaneous quantities at timestep 1
if storethQ && thTimePoints(1)
    thOut(:,:,1) = thm;
    thTimePoint = 2;
end

%store forward scheme quantities at timestep 1
if scheme == 0
    Vn = (th0-thm)/dt;
    if storeVQ && VTimePoints(1)
        VOut(:,:,1) = Vn;
        VTimePoint = 2;
    end
    if storeVtotQ && VtotTimePoints(1)
        VtotOut(:,1) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
        VtotTimePoint = 2;
    end
end


%handle central and backward scheme quantities at timestep 1
if scheme == 1 || scheme == 2
    if storeVQ && VTimePoints(1)
        VTimePoint = 2;
    end
    if storeVtotQ && VtotTimePoints(1)
        VtotTimePoint = 2;
    end
end

%handle full scheme quantities at current timestep
if storeIQ && ITimePoints(1)
    ITimePoint = 2;
end
if storeEQ && ETimePoints(1)
    ETimePoint = 2;
end

%do computations 
for i = 2:Nt-1
    
    %check if I and Temp are time dependent
    if ~TCompact(2); iTemp = i; else; iTemp = 1; end
    if ~IExtCompact(2); iI = i; else; iI = 1; end
    
    if (storeVQ && VTimePoints(i)) || (storeVtotQ && VtotTimePoints(i)) || (storeEQ && ETimePoints(i))
        computeVQ = true;
    else
        computeVQ = false;
    end
    
    %store instantaneous quantities at current timestep
    if storethQ && thTimePoints(i)
        thOut(:,:,thTimePoint) = th0;
        thTimePoint = thTimePoint + 1;
    end
    if storeEQ && ETimePoints(i)
        if cpQ
            EOut(:,ETimePoint) = shiftdim(sum(icp(Ic,th0),1),1);
        else
            EOut(:,ETimePoint) = shiftdim(sum(Ic.*(1-cos(th0)),1),1);
        end
    end
    
    %store backward scheme quantities at current timestep
    if scheme == 2
        if computeVQ
            Vn = (th0-thm)/dt;
        end
        if storeVQ && VTimePoints(i)
            VOut(:,:,VTimePoint) = Vn;
            VTimePoint = VTimePoint + 1;
        end
        if storeVtotQ && VtotTimePoints(i)
            VtotOut(:,VtotTimePoint) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
            VtotTimePoint = VtotTimePoint + 1;
        end
        if storeEQ && ETimePoints(i)
            EOut(:,ETimePoint) = EOut(:,ETimePoint) + shiftdim(sum(betaC.*Vn.^2,1),1)/2;
        end
    end

    %core computations of this timestep
    ILn = ILBase.*IExt(1,IExtprob,iI) ;
    IEn = IExtBaseJ.*IExt(1,IExtprob,iI);
    if cpQ
        yn = cp(Ic,th0);
    else
        yn = Ic.*sin(th0);
    end
    yn = yn + randn(Nj,W).*sqrt(2*T(1,Tprob,iTemp)/dt./Rn)- IEn + C0.*th0 + Cm.*thm;
    if LMode == 0
        thp = (AT*(D\(A*(g-th0+ILn)))/betaL-yn)./Cp;    
    else
        thp = (AT*(DB\(A*(g-th0+ILn)))-yn)./Cp;
    end
    
    %check if I needs to be computed
    if (storeIQ && ITimePoints(i)) || (storeEQ && ETimePoints(i))
        I = yn + IEn + Cp.*thp;
    end
    

    
    %store full scheme quantities at current timestep
    if storeIQ && ITimePoints(i)
        IOut(:,:,ITimePoint) = I;
        ITimePoint = ITimePoint + 1;
    end
    if storeEQ && ETimePoints(i)
        EOut(:,ETimePoint) = EOut(:,ETimePoint) + shiftdim(sum(I.*(betaL*I),1),1)/2;
    end
    
    %store central scheme quantities at current timestep
    if scheme == 1
        if computeVQ
            Vn = (thp-thm)/dt/2;
        end
        if storeVQ && VTimePoints(i)
            VOut(:,:,VTimePoint) = Vn;
            VTimePoint = VTimePoint + 1;
        end
        if storeVtotQ && VtotTimePoints(i)
            VtotOut(:,VtotTimePoint) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
            VtotTimePoint = VtotTimePoint + 1;
        end
        if storeEQ && ETimePoints(i)
            EOut(:,ETimePoint) = EOut(:,ETimePoint) + shiftdim(sum(betaC.*Vn.^2,1),1)/2;
        end
    end
    
    %store forward scheme quantities at current timestep
    if scheme == 0
        if computeVQ
            Vn = (thp-th0)/dt;
        end
        if storeVQ && VTimePoints(i)
            VOut(:,:,VTimePoint) = Vn;
            VTimePoint = VTimePoint + 1;
        end
        if storeVtotQ && VtotTimePoints(i)
            VtotOut(:,VtotTimePoint) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
            VtotTimePoint = VtotTimePoint + 1;
        end
        if storeEQ && ETimePoints(i)
            EOut(:,ETimePoint) = EOut(:,ETimePoint) + shiftdim(sum(betaC.*Vn.^2,1),1)/2;
        end
    end  
    
    %prepare quantities for next timestep
    if storeEQ && ETimePoints(i)
        ETimePoint = ETimePoint + 1;
    end
    thm = th0;
    th0 = thp;
end

%store instantaneous quantities at current timestep
if storethQ && thTimePoints(Nt)
    thOut(:,:,thTimePoint) = th0;
end

%store backward scheme quantities at last timestep
if scheme == 2
    Vn = (th0-thm)/dt;
    if storeVQ && VTimePoints(Nt)
        VOut(:,:,VTimePoint) = Vn;
    end
    if storeVtotQ && VtotTimePoints(Nt)
        VtotOut(:,VtotTimePoint) = sum(IExtBaseJ.*Vn,1)'/IExtBaseTotal;
    end
end
end
