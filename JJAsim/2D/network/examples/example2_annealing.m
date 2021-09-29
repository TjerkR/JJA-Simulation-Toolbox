close all; 
% - this example shows how to find low energy stationairy states with 
%   simulated annealing.
% - Done in the presence of an external magnetic field, which is represented 
%   by the frustration factor f.

%array size
Nx = 10;                 %nr of islands in x-direction
Ny = 10;                 %nr of islands in y-direction
ax = 1;                  %island spacing in x-direction
ay = 1;                  %island spacing in y-direction

%other parameters
f = 1/3;                 %frustration factor
inputMode = 'sweep';     %input mode
IExt = 0;                %External current
t = (0:0.5:20000)';      %time points
Nt = length(t);
z = 0;                   %phase zone
th1 = 0;                 %initial condition 
T = linspace(0.15,0,Nt); %annealing temperature profile

%create square array
array = JJAsim_2D_network_square(Nx,Ny);

%time evolution
out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1);
nOut = JJAsim_2D_network_method_getn(array,out.th,z);   

%visualize annealing process
selectedTimePoints = false(Nt,1);
selectedTimePoints(1:20:end) = true;
JJAsim_2D_visualize_movie(array,t,zeros(array.Np,Nt),out.I,'showPathQuantityQ',...
    true,'pathQuantity',nOut,'selectedTimePoints',selectedTimePoints,...
    'showIslandsQ',false,'arrowColor',[0,0,0],'pathQuantityLabel','n');
close all;

%visualize final configuration 
nFinal = nOut(:,:,end-1);
IFinal = out.I(:,:,end-1);
%JJAsim_2D_visualize_snapshot(array,nFinal,IFinal,'showIExtBaseQ',false);

disp(['energy of final state (per junction): ',num2str(out.E(end)/array.Nj)])






% plot(t,out.E,'LineWidth',1.5,'Color',[0,0,0]); 
% xlabel('$t$','Interpreter','latex','FontSize',15); 
% ylabel('$E(t)$','Interpreter','latex','FontSize',15); 
% title('anneal temperature profile','Interpreter','latex','FontSize',15);
% ah = gca;
% ah.LineWidth = 1.5;
% ah.FontSize = 15;
% ah.XAxis.TickLabelInterpreter = 'latex';
% ah.YAxis.TickLabelInterpreter = 'latex';
% ah.XLim = [-500,20000];
% ah.YLim = [50,170];
% 
% [n,ia,ic] = unique(squeeze(nOut)','rows');
% th1 = JJAsim_2D_network_stationairyState_approx_london(array,n',f);
% out = JJAsim_2D_network_stationairyState(array,'list',IExt,f,n',th1);
% Eu = out.E;
% Eu(~out.solutionQ) = nan;
% hold on; 
% plot(t,Eu(ic),'LineWidth',1.5,'Color',[0.8,0.2,0])
% 




