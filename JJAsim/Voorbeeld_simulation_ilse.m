close all
clear all
filename = 'test'
folder = fileparts(which(filename))
addpath(genpath(folder))
S=[];
tic
N1=100;
N2=100;
L=50;
filedirectory = '~/matlab/JJAsim-master/JJAsim-master/';
array = JJAsim_2D_network_square(N1,N2,1,1,'y');
xy = array.nodePosition;

bound1 = ((N1-1)/2)-L/2;
bound2 = ((N1-1)/2)+L/2;
bound3 = ((N2-1)/2)-L/2;
bound4 = ((N2-1)/2)+L/2;

nodeNrs = (xy(:,1) > bound1 & xy(:,1) < bound2) & (xy(:,2) > bound3 & xy(:,2) < bound4);
nodeNrs = find(nodeNrs);
I1 = array.IExtBase;
array = JJAsim_2D_network_removeNodes(array,nodeNrs);
I2 = array.IExtBase;
toc
% z = 0;
% inputMode = 'sweep';

x0 = 50;
y0 = 50;
n0 = 1;
f = 0.07;

% [th,I,n] = JJAsim_2D_network_stationairyState_approx_arctan(array,x0,y0,n0,f);
% phi = JJAsim_2D_network_method_getphi(array,th);
% 
% JJAsim_2D_visualize_snapshot(array,n,I,'showNodeQuantityQ',true,'nodeQuantity',phi);


holeNr = find(array.pathArea > 1);
Imax2=[];
hold on;
nHole1 = -50:1:50;
for kk=1:length(nHole1)
    kk
    nHole = nHole1(kk);
%nHole = -2;
nHole2 = 0;
extra = 23;


fList = -0.5:0.0001:0.5; 
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
hold off
 filename = ['N1_' num2str(N1) '_N2_' num2str(N2) '_L_' num2str(L) '.jpeg'];
saveas(fig, strcat(filedirectory,filename),'jpeg')
Imax2
  filename1 = ['N1_' num2str(N1) '_N2_' num2str(N2) '_L_' num2str(L) '_' num2str(P(1)) '_' num2str(P(2)) '.jpeg'];
 saveas(fig1, strcat(filedirectory,filename1),'jpeg')


 filename3=['~/matlab/JJAsim-master/JJAsim-master/Matrix_Iext_N1_' num2str(N1) '_N2_' num2str(N2) '_L_' num2str(L)];
save(filename3,'Imax2');
