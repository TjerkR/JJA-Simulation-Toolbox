%% Example - Calculate Ic(B) for batch of disordered arrays
% Tjerk Reintsema
% 18-10-2021
% 
% Calculates and saves Ic(B) data for a batch of given arrays.
clearvars

savedir = '/home/s1736159/MATLAB/data/';

%% Add current folder to path, start parallel pool
p = genpath(pwd);
addpath(p);

parpool

tic

%% Inputs

% Device F (only hole size, array size is actually 50 but takes too long)
N = 34;
L = 30;

variances = linspace(0, 0.10, 21);
sigmas = sqrt(variances);

Nj = 2*N*(N-1);
Ics = zeros(Nj, length(sigmas));
variances_actual = zeros(1, length(sigmas));
means_actual = zeros(1, length(sigmas));
for i = 1:length(sigmas)
    Ics(:,i) = sigmas(i) * randn(Nj,1) + 1;
    
    variances_actual(i) = var(Ics(:,i));
    means_actual(i) = mean(Ics(:,i));
end

hole_shape = 'square';

mirroring = false; % do NOT use mirroring with disordered arrays!
parQ = false;
numcores = 4;

fstart = -6E-3;
fstep = 6E-5;
fstop = 6E-3;
f_list = fstart:fstep:fstop;

%% Calculations
for j = 1:length(sigmas)
    fprintf('%da %dh %dvar %dmean\n\n',N,L,variances_actual(j),means_actual(j))
    
    nHole_list = -3:1:3; % this has to be inside the loop!!
    
    Ic = Ics(:, j);

    % Generate array geometry
    array = generate_array(N, L, 'hole_shape', hole_shape, 'Ic', Ic);

    % Calculate Ic(B) for different number of vortices in hole
    holeNr = find(array.pathArea > 1); % ??? no idea
    S = [];
    Ic_f_list = [];
    Ic_peak_index_list = []; 
    fprintf('Starting:\n\n')
    for nHole = nHole_list

        Ic_f = zeros(size(f_list));
        arrayNp = array.Np;
        parfor i = 1:length(f_list)

            f = f_list(i);

            n = zeros(arrayNp, 1); % ???
            n(holeNr) = nHole;
            Ic_out = JJAsim_2D_network_method_getIExtMax(array, n, f, ...
                'parallelQ', parQ, ...
                'cores', numcores);
            n_out = JJAsim_2D_network_method_getn(array, Ic_out.th, 0);

            if Ic_out.flag == 0
                Ic_f(i) = Ic_out.IExtLowerBound;
            end

        end

%         [Ic_peak, Ic_peak_index] = max(Ic_f);
%         Ic_peak_index_list = [Ic_peak_index_list, Ic_peak_index];
        Ic_f_list = [Ic_f_list; Ic_f];

        disp(nHole)
    end

    if mirroring
       Ic_f_list = [fliplr(flipud(Ic_f_list(2:end,:))); Ic_f_list];
       nHole_list = [-fliplr(nHole_list(2:end)), nHole_list];
    end

    Ic_f_max = max(Ic_f_list);


    %% Saving data

    filename = strcat(string(N),'a',string(L),'h_',hole_shape,'_',erase(sprintf('%.4f', variances(j)),'.'),'var');

    save(strcat(savedir, filename), 'Ic_f_max', 'Ic_f_list', 'f_list', 'nHole_list', 'Ics', 'variances_actual', 'means_actual');   

end

toc


