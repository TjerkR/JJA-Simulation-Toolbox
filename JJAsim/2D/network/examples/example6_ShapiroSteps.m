close all; 
% - magnetoresistance curves for several bias currents

%array size
Nx = 21;                
Ny = 21;

%create square array
array = JJAsim_2D_network_square(Nx,Ny);

%time
t = (0:0.1:300)';       
Nt = length(t);

%bias current components
freq = 0.25;
Amp1 = 0.5;
Amp2 = 1;
IDC = linspace(0,2,201)';

IAC1 = Amp1*sin(freq*t)';
IAC2 = Amp2*sin(freq*t)';
IExt = [IDC + IAC1;IDC + IAC2];    

%other parameters
f = 0;                  %frustration factor
inputMode = 'sweep';    %input mode (only relevant if computing multiple problems in one function call)
z = 0;                  %phase zone
th1 = 0;                %initial condition
T = 0.01;               %temperature

%time evolution (storethQ and storeIQ are set to false to reduce memory usage)
out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1,'storethQ',false,...
    'storeIQ',false);

%compute array resistance
V = mean(out.Vtot(:,round(Nt/3):end),2);
V = reshape(V,length(IDC),2);

%plot IV curve
plot(IDC,V,'LineWidth',1.5)
ah = gca;
ah.Box = 'on';
ah.LineWidth = 1.5;
xlabel('$I$','Interpreter','latex','FontSize',13)
ylabel('$V$','Interpreter','latex','FontSize',13)
legend({['Amplitude = ',num2str(Amp1)],['Amplitude = ',num2str(Amp2)]},'Location','NorthWest','Interpreter','latex','FontSize',13)
ah.FontSize = 15;
