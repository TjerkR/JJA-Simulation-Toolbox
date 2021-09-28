%% Single script to run on server


%% Adding current folder to path
p = genpath(pwd);
addpath(p);

parpool

tic

%% Inputs

L = 10;
N = 20;
hole_shape = 'square';

nHole_list = 0:1:2; 

mirroring = true; % JUST nHole_list!!!
parQ = false;
numcores = 4;

fstart = -10E-2;
fstep = 5E-4;
fstop = 10E-2;
f_list = fstart:fstep:fstop;

%% Calculations
fprintf('%da %dh\n\n',N,L)
 
% Generate array geometry
array = generate_array(N, L, hole_shape);


% Calculate Ic(B) for different number of vortices in hole

holeNr = find(array.pathArea > 1);
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

% Ic_f_max = max(Ic_f_list);

if mirroring
%    Ic_f_max = [fliplr(Ic_f_max(2:end)), Ic_f_max];
%    Ic_f_list = [fliplr(Ic_f_list(:, 2:end)), Ic_f_list];
   Ic_f_list = [fliplr(flipud(Ic_f_list(2:end,:))); Ic_f_list];
%    f_list = [-fliplr(f_list(2:end)), f_list];
   nHole_list = [-fliplr(nHole_list(2:end)), nHole_list];
end

Ic_f_max = max(Ic_f_list);


%% Saving data

savedir = '/home/s1736159/MATLAB/data/';
% savedir = '.\';

filename = strcat(string(N),'a',string(L),'h');

% save(strcat(savedir, filename), 'Ic_f_max', 'Ic_f_list', 'f_list', 'nHole_list');

toc