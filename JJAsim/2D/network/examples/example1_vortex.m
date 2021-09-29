close all; 
% - this example shows how to compute a stationairy single vortex state in  
%   a square array.
% - It does this in an array with and without an inductance parameter. 

%array size
Nx = 10;                %nr of islands in x-direction
Ny = 10;                %nr of islands in y-direction
ax = 1;                 %island spacing in x-direction
ay = 1;                 %island spacing in y-direction

%vortex configuration
x_n = 4.5;              %(x_n,y_n) is the coordinate for the placed vortex
y_n = 4.5;
n = 1;                  %n is the vorticity of the placed vortex

%other parameters
f = 0;                  %frustration factor
inputMode = 'sweep';    %input mode (only relevant if computing multiple problems)
IExtDirection = 'x';    %external current direction, either 'x' or 'y'
IExt = 0;               %External current
betaLList = [0,3];      %inductance parameter. Each junction has a self inductance 
                        % of L = betaL*Phi0/(2*pi*I0)

%compute stationairy states
for i = 1:2
    betaL = betaLList(i);
    
    %create square array
    array = JJAsim_2D_network_square(Nx,Ny,ax,ay,IExtDirection,'betaL',betaL);
    
    %make initial guess
    [th,I] = JJAsim_2D_network_stationairyState_approx_arctan(array,x_n,y_n,n,f);
    z = 0; %When using approx_arctan, the phase zone z must be zero.
    
    %compute exact stationairy state 
    out = JJAsim_2D_network_stationairyState(array,inputMode,IExt,f,z,th);
    
    %vortex configuration of stationairy state
    nOut = JJAsim_2D_network_method_getn(array,out.th,z);
    
    %compute island phases
    phiOut = mod(JJAsim_2D_network_method_getphi(array,out.th),2*pi);
    
    %display stationairy state
    figure;
    if i == 1
        title('no inductance parameter')
    else
        title(['inductance parameter $\beta_L = $',num2str(betaL)])
    end
    h = JJAsim_2D_visualize_snapshot(array,nOut,out.I,'showNodeQuantityQ',true,...
        'nodeQuantity',phiOut,'nodeColorLimits',[0,2*pi],'nodeQuantityLabel',...
        '$\phi$','showIExtBaseQ',false,'nodeDiameter',0.35,'FontSize',20,'arrowWidth',1.7);
    h.axisHandle.XAxis.Visible = 'off';
    h.axisHandle.YAxis.Visible = 'off';
end
