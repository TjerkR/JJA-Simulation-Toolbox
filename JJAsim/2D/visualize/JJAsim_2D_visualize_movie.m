function JJAsim_2D_visualize_movie(array,t,n,I,varargin)
%JJAsim_2D_movie(array,t,n,I,varargin) 
%
%DESCRIPTION
% - Displays a movie of a time evolution simulation on a 2D josephson junction array.
% - Currents I are displayed with arrows where the current is proportional to the strength
% - Vortices n are displayed as circles.
% - Currents and vortices can be displayed with a 'normal' or 'fancy' draw mode. Normal 
%   is more efficient. Set with arrowType and vortexType.
% - Optionally one can display node quantities with colored nodes. Can for example be
%   used to display node phases or node potential. 
% - Optionally one can display path quantities with colored path area. Can for example be
%   used to display path currents or self fields.
% - Optionally one can display the nodes in which external current is injected and ejected
%   with showIExtBaseQ. Injected nodes get a dot and ejection nodes get a cross.
% - Movie can be recorded into an .avi file.
%
%FIXED INPUT
% array               struct      information about Josephson junction array.
% t                   Nt by 1     time points
% n                   Np by Nt    vortex configuration
% I                   Nj by Nt    current configuration
%
%VARIABLE INPUT
% figurePosition      empty       automatic figure position
%                     1 by 1      figure is centered, between 0 and 1 where 1 = fullscreen.  
%                     1 by 4      manual figure position.      
% selectedTimePoints  Nt by 1     The entries that are true are displayed in the movie.
% FontName            string      font name
% FontSize            1 by 1      font size
% showVorticesQ       1 by 1      If true, vortices are displayed in the centre of paths.
% vortexDiameter      1 by 1      Diameter of displayed vortices.
% vortexColor         1 by 3      RGB triplet for vortex color
% antiVortexColor     1 by 3      RGB triplet for antivortex color  
% vortexType          string      'normal' or 'fancy'. Draw mode for vortices.
% showGridQ           1 by 1      draw a line grid over all junctions
% gridWidth           1 by 1      width of gridlines
% gridColor           1 by 3      RGB triplet for grid color
% showCurrentQ        1 by 1      If true, current is displayed with arrows whose length
%                                 is proportional to the magnitude of the current.
% arrowWidth          1 by 1      width scale factor of current arrays.
% arrowLength         1 by 1      length scale factor of current arrays.
% arrowColor          1 by 3      RGB triplet for arrow color  
% arrowType           string      'normal' or 'fancy'. Draw mode for current arrows.
% showIslandsQ        1 by 1      If true, nodes are shown as circles.
% nodeDiameter        1 by 1      Diameter of displayed nodes.
% nodeColor           1 by 3      RGB triplet for color of nodes 
%                                 (ignored if showIslandQuantityQ)
% showIslandQuantityQ 1 by 1      if true, nodeQuantity is displayed with colored nodes.
% nodeQuantity        Nn by Nt    quantity displayed with colored nodes.   
% nodeColorLimits     empty       automatic color limits for nodeQuantity
%                     1 by 2      manual color limits for nodeQuantity
% nodeQuantityLabel   string      Colorbar label for nodeQuantity  
% showPathQuantityQ   1 by 1      if true, pathQuantity is displayed by coloring the 
%                                 enclosed area of each path.
% pathQuantity        Np by Nt    quantity displayed with colored path areas.   
% pathColorLimits     empty       automatic color limits for pathQuantity
%                     1 by 2      manual color limits for pathQuantity
% pathQuantityLabel   string      Colorbar label for pathQuantity    
% pathQuantityAlpha   1 by 1      Transparancy of colored path areas.  
% showIExtBaseQ       1 by 1      Dispay the nodes where external current is in/ejected.  
% IExtBaseColor       1 by 3      RGB color triplet for symbols displaying external current.
% framePause          1 by 1      extra time between frames
% saveQ               1 by 1      true to generate a .avi file
% filename            string      filename of generated .avi file
% framerate           1 by 1      framerate of generated .avi file
% compression         1 by 1      compression factor of generated .avi file. Between 0 and 1.

%check if array is 2D
if array.ndims ~= 2
    error('input must be 2D array')
end

Nn = array.Nn;
Nj = array.Nj;
Np = array.Np;

t = reshape(double(t),[],1);
Nt = length(t);
n = reshape(double(n),Np,[]);
I = reshape(double(I),Nj,[]);

inputParameters = {
    'figurePosition'            [];
    'selectedTimePoints'        true(Nt,1)
    'FontName'                  'Arial'
    'FontSize'                  12
    'showVorticesQ'             true
    'vortexDiameter'            0.3
    'vortexColor'               [0,0,0]
    'antiVortexColor'           [1,0,0]
    'vortexType'                'fancy'
    'showGridQ'                 true
    'gridWidth'                 0.05
    'gridColor'                 [.8,.8,.8]
    'showCurrentQ'              true
    'arrowWidth'                1
    'arrowLength'               1
    'arrowColor'                [0,0,1]
    'arrowType'                 'fancy'
    'showIslandsQ'              true
    'nodeDiameter'            0.3
    'nodeColor'               [1,1,1]
    'showIslandQuantityQ'       false
    'nodeQuantity'            []
    'nodeColorLimits'         []
    'nodeQuantityLabel'       ''
    'showPathQuantityQ'         false
    'pathQuantity'              []
    'pathColorLimits'           []
    'pathQuantityLabel'         ''
    'pathQuantityAlpha'         1
    'showIExtBaseQ'             true
    'IExtBaseColor'             [0,0,0]
    'framePause'                0.005
    'saveQ'                     false
    'filename'                  'JJAmovie'
    'framerate'                 24
    'compression'               0.5
    };

options = JJAsim_method_parseOptions(inputParameters,varargin,'JJAsim_2D_movie');

figurePosition = options.figurePosition;
selectedTimePoints = options.selectedTimePoints;
fontSize = options.FontSize;
fontName = options.FontName;
showVorticesQ = options.showVorticesQ;
vortexDiameter = options.vortexDiameter;
vortexColor = options.vortexColor;
antiVortexColor = options.antiVortexColor;
vortexType = options.vortexType;
showGridQ = options.showGridQ;
gridWidth = options.gridWidth;
gridColor = options.gridColor;
showCurrentQ = options.showCurrentQ;
arrowWidth = options.arrowWidth;
arrowLength = options.arrowLength;
arrowColor = options.arrowColor;
arrowType = options.arrowType;
showIslandsQ = options.showIslandsQ;
nodeDiameter = options.nodeDiameter;
nodeColor = options.nodeColor;
showIslandQuantityQ = options.showIslandQuantityQ;
nodeQuantity = options.nodeQuantity;
nodeColorLimits = options.nodeColorLimits;
nodeQuantityLabel = options.nodeQuantityLabel;
showPathQuantityQ = options.showPathQuantityQ;
pathQuantity = options.pathQuantity;
pathColorLimits = options.pathColorLimits;
pathQuantityLabel = options.pathQuantityLabel;
pathQuantityAlpha = options.pathQuantityAlpha;
showIExtBaseQ = options.showIExtBaseQ;
IExtBaseColor = options.IExtBaseColor;
framePause = options.framePause;
saveQ = options.saveQ;
filename = options.filename;
framerate = options.framerate;
compression = options.compression;

if showIExtBaseQ && isfield(array,'IExtBase') == false
    error('Cannot display IExtBase because array does not contain the field IExtBase');
end
if showIslandQuantityQ && showPathQuantityQ
    error('cannot simultaneously display an node quantity and a path quantity'); 
end

%prepare input
nodePosition = reshape(array.nodePosition,[],2);
junctionPosition = reshape(array.junctionPosition,[],4);
pathCentroid = reshape(array.pathCentroid,[],2);
pathPosition = reshape(array.pathPosition,[],1);
if showIslandQuantityQ
    nodeQuantity = reshape(nodeQuantity,Nn,[]);
end
if showPathQuantityQ
    pathQuantity = reshape(pathQuantity,Np,[]);
end
if showIExtBaseQ
    IExtBase = reshape(array.IExtBase,[],1);
else
    IExtBase = [];
end
t = JJAsim_method_checkInput(t,'double',Nt,0,'t');
selectedTimePoints = JJAsim_method_checkInput(selectedTimePoints,'logical',Nt,0,'selectedTimePoints');
n = JJAsim_method_checkInput(n,'double',[Np,Nt],[0,0],'n');
I = JJAsim_method_checkInput(I,'double',[Nj,Nt],[0,0],'I');
if showIslandQuantityQ
    nodeQuantity = JJAsim_method_checkInput(nodeQuantity,'double',...
        [Nn,Nt],[0,0],'nodeQuantity');
    if ~isempty(nodeColorLimits)
        nodeColorLimits = JJAsim_method_checkInput(nodeColorLimits,'double',...
            [1,2],[1,0],'nodeColorLimits');
    end
end
if showPathQuantityQ
    pathQuantity = JJAsim_method_checkInput(pathQuantity,'double',...
        [Np,Nt],[0,0],'pathQuantity');
    if ~isempty(pathColorLimits)
        pathColorLimits = JJAsim_method_checkInput(pathColorLimits,'double',...
            [1,2],[1,0],'pathColorLimits');
    end
end

%produce movie
JJAsim_priv_movie_2D(nodePosition,junctionPosition,pathCentroid,pathPosition,t,...
    selectedTimePoints,n,...
    I,figurePosition,fontSize,fontName,showVorticesQ,vortexDiameter,vortexColor,...
    antiVortexColor,vortexType,showGridQ,gridWidth,gridColor,showCurrentQ,arrowWidth,...
    arrowLength,arrowColor,arrowType,...
    showIslandsQ,nodeDiameter,nodeColor,showIslandQuantityQ,nodeQuantity,...
    nodeColorLimits,nodeQuantityLabel,showPathQuantityQ,pathQuantity,...
    pathColorLimits,pathQuantityLabel,pathQuantityAlpha,showIExtBaseQ,IExtBase,...
    IExtBaseColor,framePause,saveQ,filename,framerate,compression);


end

