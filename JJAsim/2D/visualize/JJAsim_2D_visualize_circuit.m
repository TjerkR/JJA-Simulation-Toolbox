function JJAsim_2D_visualize_circuit(array,varargin)
%JJAsim_visualize_circuit(array,varargin)
%
%DESCRIPTION
% - Displays a circuit representation of an array with the symbols of components.
% - Inductors are not shown.
% - Nodes where external current is injected are shown with a dot and nodes where
%   external current is ejected are displayed with a cross.
%
%VARIABLE INPUT
% showQuantitiesQ   1 by 1      If true, Rn, Ic, betaC and IExtBase are visualized.
% FontSize          1 by 1      Font size
% textColor         1 by 3      RGB triplet for text color.
% lineWidth         1 by 1      Width of all lines (in pixels)
% nodeDiameter      1 by 1      Diameter of nodes (circuit nodes)
% JLength           1 by 1      Length of all junctions (as a fraction of its wire)
% JWidth            1 by 1      Width of all junctions (as a fraction of its wire)
% RWidth            1 by 1      Resistor is rectangle with this width (fraction of wire)
% RLength           1 by 1      ... and this length (fraction of wire)
% JEsize            1 by 1      Width of Josephson element as a fraction of wire
% CGap              1 by 1      Capcitor gap size as a fraction of wire
% CWidth            1 by 1      Capacitor width as a fraction of wire
    
%import optional property-value pairs
inputParameters = {
    'showQuantitiesQ'   true
    'FontSize'          18
    'textColor'         [0.05,0.5,1]
    'lineWidth'         1.5
    'nodeDiameter'      0.16
    'JLength'           0.42
    'JWidth'            0.14
    'RWidth'            0.08
    'RLength'           0.2
    'JEsize'            0.12
    'CGap'              0.04
    'CWidth'            0.14
    };
options = JJAsim_method_parseOptions(inputParameters,varargin,mfilename);

showQuantities = options.showQuantitiesQ;
FontSize = options.FontSize;
textColor = options.textColor;
lineWidth = options.lineWidth;
nodeDiameter = options.nodeDiameter;
JLength = options.JLength;
JWidth = options.JWidth;
RWidth = options.RWidth;
RLength = options.RLength;
JEsize = options.JEsize;
CGap = options.CGap;
CWidth = options.CWidth;

is = array.nodePosition;
j1 = array.junctionIsland1;
j2 = array.junctionIsland2;
Nj = array.Nj;


Rn = array.Rn;
if array.RnCompactQ
    Rn = repmat(Rn,Nj,1);
end
Ic = array.Ic;
if array.IcCompactQ
    Ic = repmat(Ic,Nj,1);
end
betaC = array.betaC;
if array.betaCCompactQ
    betaC = repmat(betaC,Nj,1);
end

elements = cell(Nj,1);
nrOfElements = zeros(Nj,1);
for i = 1:Nj
    elnr = 1;
    if Rn(i) ~= 0
        elements{i}(elnr) = 1;
        elnr = elnr + 1;
    end
    if Ic(i) ~= 0
        elements{i}(elnr) = 2;
        elnr = elnr + 1;
    end
    if betaC(i) ~= 0
        elements{i}(elnr) = 3;
        elnr = elnr + 1;
    end
    nrOfElements(i) = elnr - 1;
end

angle = linspace(0,2*pi,20);
cx = cos(angle)*nodeDiameter/2;
cy = sin(angle)*nodeDiameter/2;
cross1 = [-1,1]*nodeDiameter/4;
cross2 = [1,-1]*nodeDiameter/4;

%draw nodes
plot((cx + is(:,1))',(cy+is(:,2))','LineWidth',lineWidth,'Color',[0,0,0]);
hold on;
ind = array.IExtBase > 0;
plot((cx/3 + is(ind,1))',(cy/3 + is(ind,2))','LineWidth',lineWidth,'Color',[0,0,0]);
ind = array.IExtBase < 0;
plot((cross1 + is(ind,1))',(cross1 + is(ind,2))','LineWidth',lineWidth,'Color',[0,0,0]);
plot((cross1 + is(ind,1))',(cross2 + is(ind,2))','LineWidth',lineWidth,'Color',[0,0,0]);
if showQuantities
    for i = 1:length(ind)
        if array.IExtBase(i) ~= 0
            text(is(i,1),is(i,2),['$',num2str(abs(array.IExtBase(i))),'$'],'Interpreter','latex','FontSize',FontSize,'Color',textColor,'HorizontalAlignment','center')
        end
    end
end
L1 = zeros(2,20*Nj);
L2 = zeros(2,20*Nj);
S1 = zeros(2,2*Nj);
S2 = zeros(2,2*Nj);
S3 = zeros(2,2*Nj);
S4 = zeros(2,2*Nj);

%draw junctions
iL = 1;
iS = 1;
for i = 1:Nj
    Ne = nrOfElements(i);
    r1 = is(j1(i),:)';
    r2 = is(j2(i),:)';
    c = (r1 + r2)/2;
    L = sqrt((r1(1)-r2(1))^2 + (r1(2)-r2(2))^2);
    khat = 1/L*[r2(1)-r1(1);r2(2)-r1(2)];
    jhat = 1/L*[r1(2)-r2(2);r2(1)-r1(1)];
    r1 = r1 + nodeDiameter/2*khat;
    r2 = r2 - nodeDiameter/2*khat;
    for j = 1:Ne
        
        p = c - JLength/2*khat + ((Ne+1)/2 - j)*JWidth*jhat;
        q = c + JLength/2*khat + ((Ne+1)/2 - j)*JWidth*jhat;
        if j == 1
            P = p;
            Q = q;
        end
        switch elements{i}(j)
            case 1
                [L1(:,iL:iL+1),L2(:,iL:iL+1),S1(:,iS),S2(:,iS),S3(:,iS),S4(:,iS)] = ...
                    getResistor(p,q,RLength,RWidth,khat,jhat);
                iL = iL + 2;
                iS = iS + 1;
            case 2
                [L1(:,iL:iL+2),L2(:,iL:iL+2)] = getJE(p,q,JEsize,khat,jhat);
                iL = iL + 3;
            case 3
                [L1(:,iL:iL+3),L2(:,iL:iL+3)] = getCapacitor(p,q,CGap,CWidth,khat,jhat);
                iL = iL + 4;
        end
        
    end
    R1 = (r1+r2)/2 - JLength/2*khat;
    R2 = (r1+r2)/2 + JLength/2*khat;
    L1(:,iL) = r1;
    L2(:,iL) = R1;
    iL = iL + 1;
    L1(:,iL) = R2;
    L2(:,iL) = r2;
    iL = iL + 1;
    L1(:,iL) = P;
    L2(:,iL) = p;
    iL = iL + 1;
    L1(:,iL) = Q;
    L2(:,iL) = q;
    iL = iL + 1;
end


L1(:,iL:end) = [];
L2(:,iL:end) = [];
S1(:,iS:end) = [];
S2(:,iS:end) = [];
S3(:,iS:end) = [];
S4(:,iS:end) = [];

%draw rest of circuit
plot([L1(1,:);L2(1,:)],[L1(2,:);L2(2,:)],'LineWidth',lineWidth,'Color',[0,0,0]);
plot([S1(1,:);S2(1,:);S3(1,:);S4(1,:);S1(1,:)],[S1(2,:);S2(2,:);S3(2,:);S4(2,:);S1(2,:)],...
    'LineWidth',lineWidth,'Color',[0,0,0]);
ah = gca;
ah.Visible = 'off';
ah.DataAspectRatio = [1,1,1];
ah.PlotBoxAspectRatio = [1,1,1];
fh = gcf;
fh.Color = [1,1,1];

%display quantities
if showQuantities
    for i = 1:Nj
        Ne = nrOfElements(i);
        r1 = is(j1(i),:)';
        r2 = is(j2(i),:)';
        c = (r1 + r2)/2;
        L = sqrt((r1(1)-r2(1))^2 + (r1(2)-r2(2))^2);
        jhat = 1/L*[r1(2)-r2(2);r2(1)-r1(1)];
        for j = 1:Ne
            p = c + ((Ne+1)/2 - j)*JWidth*jhat;
            switch elements{i}(j)
                case 1
                    if length(array.Rn) == 1
                        Rdisp = array.Rn;
                    else
                        Rdisp = array.Rn(i);
                    end
                    text(p(1),p(2),['$',num2str(Rdisp),'$'],'Interpreter','latex','FontSize',FontSize,'Color',textColor,'HorizontalAlignment','center')
                case 2
                    if length(array.Ic) == 1
                        Icdisp = array.Ic;
                    else
                        Icdisp = array.Ic(i);
                    end
                    text(p(1),p(2),['$',num2str(Icdisp),'$'],'Interpreter','latex','FontSize',FontSize,'Color',textColor,'HorizontalAlignment','center')
                    
                case 3
                    if length(array.betaC) == 1
                        Cdisp = array.betaC;
                    else
                        Cdisp = array.betaC(i);
                    end
                    text(p(1),p(2),['$',num2str(Cdisp),'$'],'Interpreter','latex','FontSize',FontSize,'Color',textColor,'HorizontalAlignment','center')
            end
        end
        
    end
end  
end

function [L1,L2] = getJE(r1,r2,sz,khat,jhat)
L1 = zeros(2,3);
L2 = zeros(2,3);
c = (r1+r2)/2;
L1(:,1) = r1;
L2(:,1) = r2;
L1(:,2) = c - (khat+jhat)*sz/2;
L2(:,2) = c + (khat+jhat)*sz/2;
L1(:,3) = c - (khat-jhat)*sz/2;
L2(:,3) = c + (khat-jhat)*sz/2;
end

function [L1,L2,S1,S2,S3,S4] = getResistor(r1,r2,length,width,khat,jhat)
L1 = zeros(2,2);
L2 = zeros(2,2);
c = (r1+r2)/2;
v1 = c - khat*length/2;
v2 = c + khat*length/2;
L1(:,1) = r1;
L2(:,1) = v1;
L1(:,2) = v2;
L2(:,2) = r2;
S1 = v1 - jhat*width/2;
S2 = v2 - jhat*width/2;
S3 = v2 + jhat*width/2;
S4 = v1 + jhat*width/2;
end

function [L1,L2] = getCapacitor(r1,r2,gap,width,khat,jhat)
L1 = zeros(2,4);
L2 = zeros(2,4);
c = (r1+r2)/2;
v1 = c - khat*gap/2;
v2 = c + khat*gap/2;
L1(:,1) = r1;
L2(:,1) = v1;
L1(:,2) = v2;
L2(:,2) = r2;
L1(:,3) = v1 - width/2*jhat;
L2(:,3) = v1 + width/2*jhat;
L1(:,4) = v2 - width/2*jhat;
L2(:,4) = v2 + width/2*jhat;
end

