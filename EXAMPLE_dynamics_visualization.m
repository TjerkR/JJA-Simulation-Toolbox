%% Example - Script to investigate array dynamics
% Tjerk Reintsema
% 29-09-2021
% 
% Visualization of vortex dynamics and phases in a disordered 34x34 array
% with a 30x30 hole.
close all
clearvars

%% Inputs

N = 34;
L = 30;

% Choose square, circle or diamond
hole_shape = 'square';

x0 = N/2;
y0 = N/2;

z = 0; % INPUT?
inputMode = 'sweep';

n0 = 1;
f = 0.05;

t = (0:0.1:100)'; % INPUT
T = 0; % INPUT
IExt = 0.8; % INPUT

% Ic = 1;
Nj = 2*N*(N-1);
Ic = sqrt(0.05) * randn(Nj,1) + 1;


%% Creating array

% Generate geometry with hole
array = generate_array(N, L, 'hole_shape', hole_shape, 'Ic', Ic);

%% Calculate phases in stationary state?

% Calculations
[th, I, n] = JJAsim_2D_network_stationairyState_approx_arctan(array, ...
    x0, y0, n0, f);
phi = JJAsim_2D_network_method_getphi(array, th);

% Visualize initial state
JJAsim_2D_visualize_snapshot(array, n, I, ...
    'showNodeQuantityQ', true, ...
    'nodeQuantity', phi);

%% Calculate dynamics

% Calculations
out = JJAsim_2D_network_simulate(array, t, inputMode, IExt, T, f, z, th);

% Get vortex configuration
nsim = JJAsim_2D_network_method_getn(array, out.th, z);

% Define timepoints
stp = false(length(t), 1); % define timepoints for movie
stp(1:2:end) = true;

% Visualize dynamics
JJAsim_2D_visualize_movie(array, t, nsim(:,1,:), sin(out.th(:,1,:)), ...
    'selectedTimePoints', stp);
