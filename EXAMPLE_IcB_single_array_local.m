%% Example - Calculate Ic(B) of single array, run locally
% Tjerk Reintsema
% 29-09-2021
% 
% Calculates and optionally saves and plots Ic(B) for a single given array.
close all
clearvars

PLOT_DATA = true;
SAVE_FIGURES = true;
SAVE_DATA = true;
filestring = "10a8h";
titlestring = "I_c of a 10 by 10 array, with an 8 by 8 hole";
% xlims = 'auto';
ylims = [0, 0.25];
figno1 = 10;
figno2 = 20;

%% Add current folder to path, start parallel pool
p = genpath(pwd);
addpath(p);

% parpool

tic

%% Inputs

N = 10;
L = 8;

hole_shape = 'square';

nHole_list = 0:1:5; 

mirroring = true;
parQ = false;
numcores = 4;

fstart = -10E-2;
fstep = 5E-4;
fstop = 10E-2;
f_list = fstart:fstep:fstop;

%% Calculations
fprintf('%da %dh\n\n',N,L)
 
% Generate array geometry
array = generate_array(N, L, 'hole_shape', hole_shape);

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

%     [Ic_peak, Ic_peak_index] = max(Ic_f);
%     Ic_peak_index_list = [Ic_peak_index_list, Ic_peak_index];
    Ic_f_list = [Ic_f_list; Ic_f];

    disp(nHole)
end

if mirroring
   Ic_f_list = [fliplr(flipud(Ic_f_list(2:end,:))); Ic_f_list];
   nHole_list = [-fliplr(nHole_list(2:end)), nHole_list];
end

Ic_f_max = max(Ic_f_list);

toc

%% Plotting

if PLOT_DATA
    % All curves
    [fig1, ax1] = plot_IcB(Ic_f_list, f_list, nHole_list, 'figno', figno1, 'titlestring', titlestring, 'ylims', ylims);

    % Maximum Ic(B)
    [fig2, ax2] = plot_IcB_max(Ic_f_max, f_list, 'figno', figno2, 'titlestring', titlestring, 'ylims', ylims);

    % Saving figures
    if SAVE_FIGURES
        if not(isfolder(".\figures\"))
            mkdir(".\figures\")
        end
        exportgraphics(ax1, strcat(".\figures\", filestring, "_all.png"), "Resolution", 400)
        exportgraphics(ax2, strcat(".\figures\", filestring, "_max.png"), "Resolution", 400)
    end
end

%% Saving data

if SAVE_DATA
    if not(isfolder(".\data\"))
        mkdir(".\data\")
    end
    save(strcat(".\data\", filestring), 'Ic_f_max', 'Ic_f_list', 'f_list', 'nHole_list');
end
