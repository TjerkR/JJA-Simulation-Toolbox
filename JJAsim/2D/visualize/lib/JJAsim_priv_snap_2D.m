function [fh,ah,cb] = JJAsim_priv_snap_2D(nodePosition,junctionPosition,pathCentroid,...
    pathPosition,n,I,figurePosition,fontSize,fontName,showVorticesQ,vortexDiameter,...
    vortexColor,antiVortexColor,vortexType,showGridQ,gridWidth,gridColor,showCurrentQ,arrowWidth,...
    arrowLength,arrowColor,arrowType,...
    showIslandsQ,nodeDiameter,nodeColor,showIslandQuantityQ,nodeQuantity,...
    nodeColorLimits,nodeQuantityLabel,showPathQuantityQ,pathQuantity,...
    pathColorLimits,pathQuantityLabel,pathQuantityAlpha,showIExtBaseQ,IExtBase,...
    IExtBaseColor)

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
fh = gcf;
ah = gca;
ah.NextPlot = 'add';
cb = [];
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

%plot path quantity
if showPathQuantityQ
    cx = nan(length(pathPosition),max(cellfun(@(x) size(x,1),pathPosition)));
    cy = cx;
    for i = 1:length(pathPosition)
        cpath = pathPosition{i};
        cx(i,1:size(cpath,1)) = cpath(:,1);
        cy(i,1:size(cpath,1)) = cpath(:,2);
    end
    patch(cx',cy',pathQuantity','EdgeColor','none','FaceAlpha',pathQuantityAlpha);
end

%plot nodes
if showIslandsQ
    
    %node circles
    if showIslandQuantityQ
        patch((X + xCicle*nodeDiameter/2)',(Y + yCicle*nodeDiameter/2)',...
            nodeQuantity','EdgeColor',[0,0,0]);
    else
        patch((X + xCicle*nodeDiameter/2)',(Y + yCicle*nodeDiameter/2)',...
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
    [x,y,u,v] = JJAsim_priv_arrows_2D([X1,Y1,X2,Y2],reshape(I,[],1),nodeDiameter*1.05,arrowLength);
    switch arrowType
        case 'normal'
            quiver(x,y,u,v,'LineWidth',arrowWidth,'Color',arrowColor,'AutoScale','off');
        case 'fancy'
            p = JJAsim_priv_quiver_2D(ah,u,v,x,y,1,1,0.7*arrowWidth);
            p.FaceColor = arrowColor;
            p.EdgeColor = arrowColor;
        otherwise
            error('unrecognized arrowType')
    end
end



%display vortices
nmax = max(abs(n));
d = vortexDiameter/2;
if showVorticesQ
    switch vortexType
        case 'fancy'
            for ntype = 1:nmax
                ind = n == ntype;
                npx = px(ind);
                npy = py(ind);
                dp = linspace(0,d,2*ntype);
                if sum(ind) > 0
                    for j = ntype:-1:1
                        if j == ntype
                            patch((npx + xCicle*dp(end))',(npy + yCicle*dp(end))',vortexColor);
                        else
                            patch((npx + xCicle*dp(2*j+1))',(npy + yCicle*dp(2*j+1))',[1,1,1]);
                            patch((npx + xCicle*dp(2*j))',(npy + yCicle*dp(2*j))',vortexColor);
                        end
                    end
                end
                ind = n == -ntype;
                npx = px(ind);
                npy = py(ind);
                if sum(ind) > 0
                    for j = ntype:-1:1
                        if j == ntype
                            patch((npx + xCicle*dp(end))',(npy + yCicle*dp(end))',antiVortexColor,'EdgeColor',antiVortexColor);
                        else
                            patch((npx + xCicle*dp(2*j+1))',(npy + yCicle*dp(2*j+1))',[1,1,1]);
                            patch((npx + xCicle*dp(2*j))',(npy + yCicle*dp(2*j))',antiVortexColor,'EdgeColor',antiVortexColor);
                        end
                    end
                end
            end
        case 'normal'
            if nmax > 1
                warning('Higher order vortices (abs(n)>1) occur, for VortexType normal these are displayed as normal vortices. Use VortexType fancy to display higer order vortices.')
            end
            patch((px(n > 0) + xCicle*d)',(py(n > 0) + yCicle*d)',vortexColor);
            patch((px(n < 0) + xCicle*d)',(py(n < 0) + yCicle*d)',antiVortexColor);
    end
end

%get plot limits
xlims = [min([X;X1;X2]) - nodeDiameter,max([X;X1;X2]) + nodeDiameter];
ylims = [min([Y;Y1;Y2]) - nodeDiameter,max([Y;Y1;Y2]) + nodeDiameter];
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
            nodeColorLimits = [min(nodeQuantity),max(nodeQuantity)];
        end
        if nodeColorLimits(1) == nodeColorLimits(2)
            nodeColorLimits = nodeColorLimits + [-0.5,0.5];
        end
        ah.CLim = nodeColorLimits;
    end
    if showPathQuantityQ
        lab = pathQuantityLabel;
        if isempty(pathColorLimits)
            pathColorLimits = [min(pathQuantity),max(pathQuantity)];
        end
        if pathColorLimits(1) == pathColorLimits(2)
            pathColorLimits = pathColorLimits + [-0.5,0.5];
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
        text(pos(1)+pos(3)*1.1,ah.Position(4)+ah.Position(2)-0.05,lab,'Units','normalized','Parent', AxesH)
    end
end

%set font size and name
set(findall(fh,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(findall(fh,'-property','Interpreter'),'Interpreter','latex')
set(findall(fh,'-property','FontSize'),'FontSize',fontSize)
set(findall(fh,'-property','FontName'),'FontName',fontName)
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

