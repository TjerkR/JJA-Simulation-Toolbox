close all
clear all
N1=20;
N2=20;
L=10;
array = JJAsim_2D_network_square(N1,N2,1,1,'y');
xy = array.nodePosition;

array = JJAsim_2D_network_square(N1,N2,1,1,'y');
xy = array.nodePosition;

% bound1 = ((N2-1)/3)-L/2;
% bound2 = ((N2-1)/3)+L/2;
% bound3 = ((N2-1)/3)-L/2;
% bound4 = ((N2-1)/3)+L/2;

bound1 = ((N1-1)/2)-L/2;
bound2 = ((N1-1)/2)+L/2;
bound3 = ((N2-1)/2)-L/2;
bound4 = ((N2-1)/2)+L/2;

% bound1 =  ((N1-1)/2)-L/2;
% bound2 = ((N1-1)/2)+L/2;
% bound3 = ((5*(N2-1)/7)-8*L/18)+1;
% bound4 = ((5*(N2-1)/7)-8*L/18)+2;
% bound1 = 1;
% bound2 = 6;
% bound3= 2;
% bound4 = 7;

nodeNrs = (xy(:,1) > bound1 & xy(:,1) < bound2) & (xy(:,2) > bound3 & xy(:,2) < bound4);
nodeNrs = find(nodeNrs);
I1 = array.IExtBase;
array = JJAsim_2D_network_removeNodes(array,nodeNrs);
I2 = array.IExtBase;

z = 0;
inputMode = 'sweep';

x0 = N1/2;
y0 = 2*N2/3;
n0 = 1;
f = 0.07;

[th,I,n] = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f);
phi = JJAsim_2D_network_method_getphi(array,th);

JJAsim_2D_visualize_snapshot(array,n,I,'showNodeQuantityQ',true,'nodeQuantity',phi);
% 
t = (0:0.1:500)';
T=0;
IExt=0.25;
th1=th;
out2 = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1);

nsim = JJAsim_2D_network_method_getn(array,out2.th,z);
stp = false(length(t),1);
stp(1:2:end) = true;

JJAsim_2D_visualize_movie(array,t,nsim(:,1,:),sin(out2.th(:,1,:)),'selectedTimePoints',stp);
JJAsim_2D_network_method_getIExtMax(array,n,f);



%% Plot the energy over time, taking account for the hole
f = linspace(0,2,101);
IExt = [0 0.2 0.4 0.6 0.8 1 1.2]';
T = [0 0.2 0.4 0.6 0.8];
partitions = 7;
for iii=1:length(T)
    T2=T(iii);
    magnetoresistance(array,t,inputMode,IExt,T2,f,z,th1,L,partitions,N1)
end
% 
%  
