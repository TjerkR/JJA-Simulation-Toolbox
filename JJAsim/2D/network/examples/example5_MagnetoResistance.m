 % - magnetoresistance curves for several bias currents

%array size
Nx = 20;               
Ny = 20;

%create square array
array = JJAsim_2D_network_square(Nx,Ny);

%other parameters
f = linspace(0,2,101);  %frustration factor
inputMode = 'sweep';    %input mode 
IExt = [0.1,0.3,0.5,0.8,1.2]';  %Bias current
t = (0:0.1:500)';       %time points
Nt = length(t);
z = 0;                  %phase zone
th1 = 0;                %initial condition
T = 0.1;                %temperature
parallelQ = true;       %if true, time evolution is computed on multiple cores.
partitions  = 4;        

%time evolution (storethQ and storeIQ are set to false to reduce memory)
out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1,'storethQ',false,...
    'storeIQ',false,'storeVQ',false,'parallelQ',parallelQ,'computePartitions',partitions);

%compute array resistance
V = mean(out.Vtot(:,round(Nt/3):end),2)/(Nx-1);
V = reshape(V,length(IExt),length(f));
R = V./IExt;

%plot RT curve
plot(f,R,'LineWidth',1.5)
xlabel('$f$','Interpreter','latex','FontSize',12)
ylabel('$R$','Interpreter','latex','FontSize',12)
lgd = cell(1,length(IExt));
for i = 1:length(IExt)
    lgd{i} = ['$I = ',num2str(IExt(i)),'$'];
end
legend(lgd,'Location','NorthEastOutside','Interpreter','latex','FontSize',12);


