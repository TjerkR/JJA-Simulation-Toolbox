%% Example - Generate and visualize disordered Ic configuration
% Tjerk Reintsema
% 27-10-2021
%
% Generates and visualizes a specific disordered Ic configuration.
% close all
clearvars
rng shuffle

%% Array and disorder setup
N = 28;
L = 20;

l_disorder = 4;

variance = 0.10;
sigma = sqrt(variance);
mu = 1;

variance_in = 0.001;
sigma_in = sqrt(variance_in);

Nj = 2*N*(N-1);

%% Generating
[Ic, array, sigma_in_actual, mu_in, mu_in_actual] ... 
    = generate_disordered_Ic(N, sigma, ...
                                "mu", mu, ...
                                "l_disorder", l_disorder, ...
                                "sigma_in", sigma_in);


%% Visualization
visualize_Ic(Ic, "L", L, "figno", 10);
visualize_simple_Ic(array, "L", L, "figno", 20);