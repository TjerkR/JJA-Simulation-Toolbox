close all
clear all
N1=30;
N2=30;
L=18;
array = JJAsim_2D_network_square(N1,N2,1,1,'y');
xy = array.nodePosition;
% filedirectory = 'C:\Users\s1530879\Documents\MATLAB\JJAsim-master\JJAsim-master';
% array = JJAsim_2D_network_square(N1,N2,1,1,'y');
% xy = array.nodePosition;

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
y0 = N1/2;
n0 = 1;
f = 0.07;

% [th,I,n] = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f);
% phi = JJAsim_2D_network_method_getphi(array,th);
% hf3 = figure(1)
% JJAsim_2D_visualize_snapshot(array,n,I,'showNodeQuantityQ',true,'nodeQuantity',phi);
% saveas(hf3,'Snapshot of sample TI dots 2','epsc')

holeNr = find(array.pathArea > 1);
S=[];
Imax2=[];
hold on;
nHole1 = -25:5:25;
for kk=1:length(nHole1)
    kk
    nHole = nHole1(kk);

fList = -0.06:0.005:0.06; 
Imax = zeros(size(fList));
for i = 1:length(fList)
    f =fList(i);
    %f
    n = zeros(array.Np,1);
    n(holeNr) = nHole;
    %n(extra) = nHole2;
    out = JJAsim_2D_network_method_getIExtMax(array,n,f);
    nOut = JJAsim_2D_network_method_getn(array,out.th,0);
    %JJAsim_2D_visualize_snapshot(array,nOut,sin(out.th))
    if out.flag == 0
        Imax(i) = out.IExtLowerBound;
    end
    
end
[s1, s2] = max(Imax);
Imax2 = [Imax2;Imax];
%fList
%Imax
hold on
fig = figure(1);
plot(fList,Imax)
title('Maximale IExt')
xlabel('f')
ylabel('Max IExt')


 a=[nHole; fList(s2)];
 S=[S a];
 P = polyfit(S(1,:),S(2,:),1);
 yfit = P(1)*S(1,:)+P(2);
 fig1 = figure(2);
 plot(S(1,:),S(2,:),'o',S(1,:),yfit)
 title('Top of IExt vs number vortices in hole')
 xlabel('Number of vortices in the hole')
 ylabel('f for which IExt is max')
end
% hold off
% filename = ['Simulation_sample_B_Au_dots_Iext_extra.jpeg'];
% saveas(fig, strcat(filedirectory,filename),'jpeg')
% Imax2
% filename1 = ['Simulation_sample__B_Au_dots_level_max_' num2str(P(1)) '_' num2str(P(2)) '_extra.jpeg'];
% saveas(fig1, strcat(filedirectory,filename1),'jpeg')
% filename2 = 'C:\Users\s1530879\Documents\MATLAB\JJAsim-master\JJAsim-master\Device_TI_dots_2\Imax_sample_A_TI_dots_2_extra.mat';
% save(filename,'Imax2');
% 
% filename3=['Matrix_Iext_sample_B_Au_dots_extra'];
% save(filename3,'Imax2');