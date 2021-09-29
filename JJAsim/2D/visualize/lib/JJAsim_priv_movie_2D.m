function JJAsim_priv_movie_2D(nodePosition,junctionPosition,pathCentroid,pathPosition,...
    t,selectedTimePoints,n,I,figurePosition,fontSize,fontName,showVorticesQ,vortexDiameter,...
    vortexColor,antiVortexColor,vortexType,showGridQ,gridWidth,gridColor,showCurrentQ,...
    arrowWidth,arrowLength,arrowColor,...
    arrowType,showIslandsQ,nodeDiameter,nodeColor,showIslandQuantityQ,nodeQuantity,...
    nodeColorLimits,nodeQuantityLabel,showPathQuantityQ,pathQuantity,pathColorLimits,...
    pathQuantityLabel,pathQuantityAlpha,showIExtBaseQ,IExtBase,IExtBaseColor,framePause,...
    saveQ,filename,framerate,compression)

%get coordinates
xCicle = cos(linspace(0,2*pi,20));
yCicle = sin(linspace(0,2*pi,20));
X = nodePosition(:,1);
Y = nodePosition(:,2);
X1 = junctionPosition(:,1);
Y1 = junctionPosition(:,2);
X2 = junctionPosition(:,3);
Y2 = junctionPosition(:,4);
px = pathCentroid(:,1);
py = pathCentroid(:,2);
normXY = sqrt((X2-X1).^2 + (Y2-Y1).^2);

%prepare figure
fh = figure;
ah = axes(fh);
hold on;
if isempty(figurePosition) || length(figurePosition) == 1
    if length(figurePosition) == 1
        scaleFactor = figurePosition;
    else
        scaleFactor = 0.5;
    end
    screensize = get(groot,'Screensize');
    W = round(screensize(3:4)*scaleFactor);
    start = screensize(1:2) + round(screensize(3:4)/2 -W/2);
    fh.Position = [start,W];
else
    fh.Position = figurePosition;
end

%plot grid
if showGridQ
    unitX = (X2-X1)./normXY;
    unitY = (Y2-Y1)./normXY;
    XGrid = [X1,X2,X2,X1] - gridWidth/2*[1 1 -1 -1].*unitY;
    YGrid = [Y1,Y2,Y2,Y1] + gridWidth/2*[1 1 -1 -1].*unitX;
    patch(XGrid',YGrid',gridColor,'EdgeColor','none')
end


%plot path quantity first timestep
if showPathQuantityQ 
    cx = nan(length(pathPosition),max(cellfun(@(x) size(x,1),pathPosition)));
    cy = cx;
    for i = 1:length(pathPosition)
        cpath = pathPosition{i};
        cx(i,1:size(cpath,1)) = cpath(:,1);
        cy(i,1:size(cpath,1)) = cpath(:,2);
    end
    pathQuantHandle = patch(cx',cy',pathQuantity(:,1)','EdgeColor','none',...
        'FaceAlpha',pathQuantityAlpha);
end

%plot nodes first timestep
if showIslandsQ
    
    %node circles
    if showIslandQuantityQ
        isHandle = patch((X + xCicle*nodeDiameter/2)',(Y + yCicle*nodeDiameter/2)',...
            nodeQuantity(:,1)','EdgeColor',[0,0,0]);
    else
        isHandle = patch((X + xCicle*nodeDiameter/2)',(Y + yCicle*nodeDiameter/2)',...
            nodeColor,'EdgeColor',[0,0,0]);
    end
    
    %plot IExtBase symbols
    if showIExtBaseQ
        ind = IExtBase < 0;
        plot((X(ind)+[-1,1]*nodeDiameter/4)',(Y(ind)+[-1,1]*nodeDiameter/4)',...
            'Color',IExtBaseColor,'LineWidth',1);
        plot((X(ind)+[-1,1]*nodeDiameter/4)',(Y(ind)+[1,-1]*nodeDiameter/4)',...
            'Color',IExtBaseColor,'LineWidth',1);
        ind = IExtBase > 0;
        patch((X(ind) + xCicle*nodeDiameter/6)',(Y(ind) + yCicle*nodeDiameter/6)',...
            IExtBaseColor,'EdgeColor',IExtBaseColor);
    end
end

%plot junctions
if showCurrentQ
    [x,y,u,v] = JJAsim_priv_arrows_2D([X1,Y1,X2,Y2],reshape(I(:,1),[],1),nodeDiameter*1.05,arrowLength);
    switch arrowType
        case 'normal'
            juncHandle = quiver(x,y,u,v,'LineWidth',arrowWidth,'Color',arrowColor,'AutoScale','off');
        case 'fancy'
            juncHandle = JJAsim_priv_quiver_2D(ah,u,v,x,y,1,1,0.7*arrowWidth);
            juncHandle.FaceColor = arrowColor;
            juncHandle.EdgeColor = arrowColor;
        otherwise
            error('unrecognized arrowType')
    end
end

% %display vortices at first timestep
d = vortexDiameter/2;
if showVorticesQ
    switch vortexType
        case 'fancy'
            [nslots,nStart,nRemoveInds,nAddValues] = getndata(n(:,selectedTimePoints));
            Ns = length(nslots);
            npatches = cell(1,Ns);
            for i = 1:Ns
                ntype = nslots(i);
                npatches{i} = gobjects(2*abs(ntype)-1,1);
                ind = nStart{i};
                npx = [-100;px(ind)];
                npy = [-100;py(ind)];
                dp = linspace(0,d,2*abs(ntype));
                if ntype > 0
                    color = vortexColor;
                else
                    color = antiVortexColor;
                end
                for j = abs(ntype):-1:1
                    if j == abs(ntype)
                        npatches{i}(2*j-1) = patch((npx + xCicle*dp(end))',(npy + yCicle*dp(end))',color);
                        npatches{i}(2*j-1).XData(:,1) =  npatches{i}(2*j-1).XData(:,1)*0;
                        npatches{i}(2*j-1).YData(:,1) =  npatches{i}(2*j-1).YData(:,1)*0;
                    else
                        npatches{i}(2*j-1) = patch((npx + xCicle*dp(2*j+1))',(npy + yCicle*dp(2*j+1))',[1,1,1]);
                        npatches{i}(2*j) = patch((npx + xCicle*dp(2*j))',(npy + yCicle*dp(2*j))',color);
                        npatches{i}(2*j-1).XData(:,1) =  npatches{i}(2*j-1).XData(:,1)*0;
                        npatches{i}(2*j-1).YData(:,1) =  npatches{i}(2*j-1).YData(:,1)*0;
                        npatches{i}(2*j).XData(:,1) =  npatches{i}(2*j).XData(:,1)*0;
                        npatches{i}(2*j).YData(:,1) =  npatches{i}(2*j).YData(:,1)*0;
                    end
                end
            end
            
        case 'normal'
            nmax = max(max(abs(n)));
            if nmax > 1
                warning('Higher order vortices (abs(n)>1) occur, for VortexType normal these are displayed as normal vortices. Use VortexType fancy to display higer order vortices.')
            end
            [nRemoveIndsVort,nRemoveIndsAnti,nAddValuesVort,nAddValuesAnti] = getndata_simple(n);
            npatchv = patch(([-100;px(n(:,1) > 0)] + xCicle*d)',([-100;py(n(:,1) > 0)] + yCicle*d)',vortexColor);
            npatcha = patch(([-100;px(n(:,1) < 0)] + xCicle*d)',([-100;py(n(:,1) < 0)] + yCicle*d)',antiVortexColor);
            npatchv.XData(:,1) = npatchv.XData(:,1)*0;
            npatchv.YData(:,1) = npatchv.YData(:,1)*0;
            npatcha.XData(:,1) = npatchv.XData(:,1)*0;
            npatcha.YData(:,1) = npatchv.YData(:,1)*0;
    end
end

%get plot limits
xlims = [min(min([X;X1;X2])) - nodeDiameter,max(max([X;X1;X2])) + nodeDiameter];
ylims = [min(min([Y;Y1;Y2])) - nodeDiameter,max(max([Y;Y1;Y2])) + nodeDiameter];
xlim(xlims)
ylim(ylims)
ah.DataAspectRatio = [1,1,1];

%set colorbar and map
if showIslandQuantityQ || showPathQuantityQ    
    
    %define colorbar and colormap
    cb = colorbar;
    colormap(parula(100));
    
    %define color limits
    if showIslandQuantityQ
        lab = nodeQuantityLabel;
    if isempty(nodeColorLimits)
        nodeColorLimits = [min(min(nodeQuantity)),max(max(nodeQuantity))];
    end
    ah.CLim = nodeColorLimits;
    end
    if showPathQuantityQ
        lab = pathQuantityLabel;
        if isempty(pathColorLimits)
            pathColorLimits = [min(min(pathQuantity)),max(max(pathQuantity))];
        end
        ah.CLim = pathColorLimits;
    end

    pos = plotboxpos(ah);
    cb.Position(2) = 0.5;
    cb.Position(1) = pos(1)+pos(3)*1.1;
    
    if isempty(lab)
        cb.Position(4) = ah.Position(4)+ah.Position(2)-0.5;
    else
        cb.Position(4) = ah.Position(4)+ah.Position(2)-0.6;
        AxesH = axes('Units', 'normalized', 'Position', [0,0,1,1], 'visible', 'off', ...
            'YLimMode', 'manual', 'YLim',  [0, 1], ...
            'XTick',    [],       'YTick', [], ...
            'NextPlot', 'add', ...
            'HitTest',  'off');
        text(pos(1)+pos(3)*1.1,ah.Position(4)+ah.Position(2)-0.05,lab,'Units','normalized','Parent', AxesH,'Interpreter','latex')
    end
end

%set title, fonsize and name
th = title(['{\itt} = ',num2str(t(1),'%0.3f')],'FontWeight','normal');
set(findall(fh,'-property','FontSize'),'FontSize',fontSize)
set(findall(fh,'-property','FontName'),'FontName',fontName)

%initialize avi file if saveQ is enabled
if saveQ
    vidObj = VideoWriter([filename,'.avi']);
    vidObj.FrameRate = framerate;
    vidObj.Quality = round(compression*100);
    open(vidObj);
end

cumPoints = cumsum(selectedTimePoints);

%start movie
for i = 2:length(t)
    ip = cumPoints(i);
    
    if selectedTimePoints(i)
        %update time value in title
        th.String = ['{\itt} = ',num2str(t(i),'%0.3f')];
        
        %update path quantity
        if showPathQuantityQ
            pathQuantHandle.CData = pathQuantity(:,i);
        end
        
        %update node quantity
        if showIslandQuantityQ
            isHandle.CData = nodeQuantity(:,i);
        end
        
        %update currents
        if showCurrentQ
            [x,y,u,v] = JJAsim_priv_arrows_2D([X1,Y1,X2,Y2],I(:,i),nodeDiameter*1.05,arrowLength);
            switch arrowType
                case 'normal'
                    juncHandle.XData = x;
                    juncHandle.YData = y;
                    juncHandle.UData = u;
                    juncHandle.VData = v;
                case 'fancy'
                    [~,x,y] = JJAsim_priv_quiver_2D(ah,u,v,x,y,1,1,0.7*arrowWidth,false);
                    juncHandle.XData = x;
                    juncHandle.YData = y;
                    
                otherwise
                    error('unrecognized arrowType')
            end
        end
        
        if showVorticesQ
            switch vortexType
                case 'fancy'
                    for k = 1:Ns
                        ntype = nslots(k);
                        npatch = npatches{k};
                        nRemove = nRemoveInds{k,ip-1};
                        nAdd = nAddValues{k,ip-1};
                        npx = px(nAdd);
                        npy = py(nAdd);
                        
                        %removing vortices
                        if ~isempty(nRemove)
                            for j = 1:2*abs(ntype)-1
                                deletePatchData(npatch(j),nRemove+1)
                            end
                        end
                        
                        %adding vortices
                        if ~isempty(nAdd)
                            dp = linspace(0,d,2*abs(ntype));
                            for j = abs(ntype):-1:1
                                if j == abs(ntype)
                                    addPatchData(npatch(2*j-1),(npx + xCicle*dp(end))',(npy + yCicle*dp(end))');
                                else
                                    addPatchData(npatch(2*j-1),(npx + xCicle*dp(2*j+1))',(npy + yCicle*dp(2*j+1))');
                                    addPatchData(npatch(2*j),(npx + xCicle*dp(2*j))',(npy + yCicle*dp(2*j))');
                                end
                            end
                        end
                        
                    end
  
            case 'normal'
                deletePatchData(npatchv,nRemoveIndsVort{i-1}+1)
                deletePatchData(npatcha,nRemoveIndsAnti{i-1}+1)
                addPatchData(npatchv,(px(nAddValuesVort{i-1}) + xCicle*d)',(py(nAddValuesVort{i-1}) + yCicle*d)');
                addPatchData(npatcha,(px(nAddValuesAnti{i-1}) + xCicle*d)',(py(nAddValuesAnti{i-1}) + yCicle*d)');
            end
        end
        drawnow;
        
        %save current frame to .avi file if saveQ is enabled
        if saveQ
            cframe = getframe(gcf);
            %cframe.cdata= imresize(cframe.cdata,0.5);
            writeVideo(vidObj,cframe);
        end
        pause(framePause);
        
    end
end
if saveQ
    close(vidObj);
end
end

function addPatchData(p,x,y)
sz = size(p.XData);
dx = size(x,2);
p.XData = [p.XData,x];
if size(p.YData,2) == sz(2)
    p.YData = [p.YData,y];
elseif size(p.YData,2) == sz(2) + dx
    p.YData(:,end-dx+1:end) = y;
else
    error('inconsistent patch data size')
end
end

function deletePatchData(p,ind)
X = p.XData;
Y = p.YData;
X(:,ind) = [];
Y(:,ind) = [];
p.XData = X;
p.YData = Y;
end

function [nslots,nStart,nRemoveInds,nAddValues] = getndata(n)
%input
%n  Np by Nt
%output
%nslots         Ns by 1          list of integers corresponding to the existing vortex types
%nStart         Ns by 1 cell     each cell is an array of all vortices of that slot in first timestep
%nRemoveInds    Ns by Nt-1 cell  which indices must be removed in advancing timestep
%nAddValues     Ns by Nt-1 cell  which vortices must be added when advancing timestep

%figure out which vortex types (slots) occur in n.
nmax = max(max(abs(n)));
occurs = false(2*nmax,1);
p = 1;
npossibilities = [-nmax:-1,1:nmax];
for ntry = [-nmax:-1,1:nmax]
    occurs(p) = sum(sum(n==ntry));
    p = p + 1;
end
nslots = npossibilities(occurs);

%initialize variables
Ns = length(nslots);
Nt = size(n,2);
nStart = cell(Ns,1);
nRemoveInds = cell(Ns,Nt-1);
nAddValues = cell(Ns,Nt-1);
nAll = cell(Ns,Nt);

%get nStart, nAll and nNrOfElements
for i = 1:Ns
   nStart{i} = find(n(:,1) == nslots(i));
   for j = 1:Nt
       nAll{i,j} = find(n(:,j) == nslots(i));
   end
end

%get nRemoveInds and nAddValues
nSim = nStart;
for i = 1:Ns
    for j = 1:Nt-1
        prev = nSim{i};
        next = nAll{i,j+1};
        [~,nRemoveInds{i,j}] = setdiff(prev,next);
        nAddValues{i,j} = setdiff(next,prev);
        prev(nRemoveInds{i,j}) = [];
        prev = [prev;nAddValues{i,j}]; %#ok<AGROW>
        nSim{i} = prev;
    end
end
end


function [nRemoveIndsVort,nRemoveIndsAnti,nAddValuesVort,nAddValuesAnti] = getndata_simple(n)
%input
%n  Np by Nt
%output
%nRemoveIndsVort    Nt-1 by 1 cell  which vortex indices must be removed in advancing timestep
%nRemoveIndsAnti    Nt-1 by 1 cell  which antivortex indices must be removed in advancing timestep
%nAddValuesVort     Nt-1 by 1 cell  which vortices must be added when advancing timestep
%nAddValuesAnti     Nt-1 by 1 cell  which antivortices must be added when advancing timestep

%initialize variables
Nt = size(n,2);
nRemoveIndsVort = cell(Nt-1,1);
nRemoveIndsAnti = cell(Nt-1,1);
nAddValuesVort = cell(Nt-1,1);
nAddValuesAnti = cell(Nt-1,1);

%get nRemoveInds and nAddValues
prevv = find(n(:,1)>0);
preva = find(n(:,1)<0);
for j = 1:Nt-1
    nextv = find(n(:,j+1)>0);
    nexta = find(n(:,j+1)<0);
    [~,nRemoveIndsVort{j}] = setdiff(prevv,nextv);
    [~,nRemoveIndsAnti{j}] = setdiff(preva,nexta);
    nAddValuesVort{j} = setdiff(nextv,prevv);
    nAddValuesAnti{j} = setdiff(nexta,preva);
    prevv(nRemoveIndsVort{j}) = [];
    preva(nRemoveIndsAnti{j}) = [];
    prevv = [prevv;nAddValuesVort{j}]; %#ok<AGROW>
    preva = [preva;nAddValuesAnti{j}]; %#ok<AGROW>
end
end



function pos = plotboxpos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%   pos:    four-element position vector, in same units as h
%
% Copyright 2010 Kelly Kearney

% Check input
if nargin < 1
    h = gca;
end
if ~ishandle(h) || ~strcmp(get(h,'type'), 'axes')
    error('Input must be an axis handle');
end

% Get position of axis in pixels
currunit = get(h, 'units');
set(h, 'units', 'pixels');
axisPos = get(h, 'Position');
set(h, 'Units', currunit);

% Calculate box position based axis limits and aspect ratios
darismanual  = strcmpi(get(h, 'DataAspectRatioMode'),    'manual');
pbarismanual = strcmpi(get(h, 'PlotBoxAspectRatioMode'), 'manual');

if ~darismanual && ~pbarismanual
    pos = axisPos;
else
    xlim = get(h, 'XLim');
    ylim = get(h, 'YLim');
    
    % Deal with axis limits auto-set via Inf/-Inf use
    if any(isinf([xlim ylim]))
        hc = get(h, 'Children');
        hc(~arrayfun( @(h) isprop(h, 'XData' ) & isprop(h, 'YData' ), hc)) = [];
        xdata = get(hc, 'XData');
        if iscell(xdata)
            xdata = cellfun(@(x) x(:), xdata, 'uni', 0);
            xdata = cat(1, xdata{:});
        end
        ydata = get(hc, 'YData');
        if iscell(ydata)
            ydata = cellfun(@(x) x(:), ydata, 'uni', 0);
            ydata = cat(1, ydata{:});
        end
        isplotted = ~isinf(xdata) & ~isnan(xdata) & ...
                    ~isinf(ydata) & ~isnan(ydata);
        xdata = xdata(isplotted);
        ydata = ydata(isplotted);
        if isempty(xdata)
            xdata = [0 1];
        end
        if isempty(ydata)
            ydata = [0 1];
        end
        if isinf(xlim(1))
            xlim(1) = min(xdata);
        end
        if isinf(xlim(2))
            xlim(2) = max(xdata);
        end
        if isinf(ylim(1))
            ylim(1) = min(ydata);
        end
        if isinf(ylim(2))
            ylim(2) = max(ydata);
        end
    end

    dx = diff(xlim);
    dy = diff(ylim);
    dar = get(h, 'DataAspectRatio');
    pbar = get(h, 'PlotBoxAspectRatio');

    limDarRatio = (dx/dar(1))/(dy/dar(2));
    pbarRatio = pbar(1)/pbar(2);
    axisRatio = axisPos(3)/axisPos(4);

    if darismanual
        if limDarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/limDarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * limDarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    elseif pbarismanual
        if pbarRatio > axisRatio
            pos(1) = axisPos(1);
            pos(3) = axisPos(3);
            pos(4) = axisPos(3)/pbarRatio;
            pos(2) = (axisPos(4) - pos(4))/2 + axisPos(2);
        else
            pos(2) = axisPos(2);
            pos(4) = axisPos(4);
            pos(3) = axisPos(4) * pbarRatio;
            pos(1) = (axisPos(3) - pos(3))/2 + axisPos(1);
        end
    end
end

% Convert plot box position to the units used by the axis
hparent = get(h, 'parent');
hfig = ancestor(hparent, 'figure'); % in case in panel or similar
currax = get(hfig, 'currentaxes');
temp = axes('Units', 'Pixels', 'Position', pos, 'Visible', 'off', 'parent', hparent);
set(temp, 'Units', currunit);
pos = get(temp, 'position');
delete(temp);
set(hfig, 'currentaxes', currax);
end


