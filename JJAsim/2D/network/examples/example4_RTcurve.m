% - Computes an RT curve (resistance versus temperature)

%array size
Nx = 20;                
Ny = 20;

%create square array
array = JJAsim_2D_network_square(Nx,Ny);

%other parameters
f = 0;                  %frustration factor
inputMode = 'sweep';    %input mode
IExt = 0.1;            %Bias current
t = (0:0.1:500)';      %time points
Nt = length(t);
z = 0;                  %phase zone
th1 = 0;                %initial condition
T = linspace(0,2,100)';   %temperature list

%time evolution (storethQ and storeIQ are set to false to reduce memory)
out = JJAsim_2D_network_simulate(array,t,inputMode,IExt,T,f,z,th1,...
    'storethQ',false,'storeIQ',false,'storeVQ',false);

%compute array resistance
R = mean(out.Vtot(:,round(Nt/3):end),2)/IExt/(Nx-1);

%plot RT curve
plot(T,R)
xlabel('T')
ylabel('R')
