close all; 
% - this example shows the time evolution on a current-biassed disordered 
%   square array where some islands are removed. For reference it also
%   shows a biassed normal square array.

%array size
Nx = 16;                  %array size
Ny = 16;
nodeRemoveFraction = 0.2; %fraction of nodes to remove

%other parameters
f = 0.2;                  %frustration factor
inputMode = 'sweep';      %input mode 
IExt = 0.5;               %External current
t = (0:0.1:100)';         %time points
Nt = length(t);
th1 = 0;                  %initial condition
z = 0;                    %phase zone
T = 0;                    %temperature

%create disordered square array
array = JJAsim_2D_network_square(Nx,Ny);
removeNodes = find(rand(array.Nn,1) < nodeRemoveFraction);
array = JJAsim_2D_network_removeNodes(array,removeNodes);

%time evolution
out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1);

%compute output vortex configuration
nOut = JJAsim_2D_network_method_getn(array,out.th,z);   

%visualize annealing process
selectedTimePoints = false(Nt,1);
selectedTimePoints(1:2:end) = true;
JJAsim_2D_visualize_movie(array,t,nOut,out.I,'selectedTimePoints',...
    selectedTimePoints,'arrowWidth',0.6);

%visualize final configuration (note, this is timestep Nt-1 because at 
%     timestep Nt one cannot obtain the current because it depends on 
%     the voltage derivative which is forward difference)
figure; 
nFinal = nOut(:,:,end-1);
IFinal = out.I(:,:,end-1);
JJAsim_2D_visualize_snapshot(array,nFinal,IFinal,'arrowWidth',0.6);
