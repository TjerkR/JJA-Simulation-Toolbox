%% Script to investigate array dynamics
% Tjerk Reintsema
% 22-03-2021
close all
clearvars

%% Inputs

N = 10;
L = 2;

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


Nj = 2*N*(N-1);
variation = 0.1 * (rand(Nj,1)-0.5);
Ic = 1*(ones(Nj,1) + variation);

%% Creating array

% Generate geometry with hole
array = generate_array(N, L, hole_shape, 'Ic', Ic);

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
