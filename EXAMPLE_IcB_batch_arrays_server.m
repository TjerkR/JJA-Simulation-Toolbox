%% Example - Calculate Ic(B) for batch of arrays, run on server
% Tjerk Reintsema
% 30-09-2021
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

Ls = (50:5:80)';
Ns = Ls+5;

A = [Ns, Ls];

hole_shape = 'square';

% nHole_list = 0:1:5; % Has to be inside the loop!

mirroring = true;
parQ = false;
numcores = 4;

fstart = 0;
fstep = 1E-5;
fstop = 1E-3;
f_list = fstart:fstep:fstop;

%% Calculations
for j = 1:length(A)
    N = A(j,1);
    L = A(j,2);
    fprintf('%da %dh\n\n',N,L)
    
    nHole_list = 0:1; % this has to be inside the loop!!

    % Generate array geometry
    array = generate_array(N, L, hole_shape);

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

    filename = strcat(string(N),'a',string(L),'h_',hole_shape);

    save(strcat(savedir, filename), 'Ic_f_max', 'Ic_f_list', 'f_list', 'nHole_list');   

end

toc


