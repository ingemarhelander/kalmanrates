
% test_poly_fit_ml - Ritar ML-polynom av olika grad för slumpad testdata.
% Kör från projektroten med:  test_poly_fit_ml
clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

%% Generate random mock data

m = 15;

% Randomly distributed step lengths
alpha = sort(0.02 + 0.48*rand(m,1));

% True underlying polynomial
F_true = 1 + 2*alpha - 1.5*alpha.^2 + 0.8*alpha.^3;

% True derivative
S_true = 2 - 3*alpha + 2.4*alpha.^2;

% Add random measurement noise
sigmaF = 0.06;
sigmaG = 0.15;

F = F_true + sigmaF*randn(m,1);
S = S_true + sigmaG*randn(m,1);

%% Test mult_q_poly_fits

% Maximum polynomial degree
qMax = 5;

% Run comparison
fits = mult_q_poly_fits(alpha, F, S, qMax);
