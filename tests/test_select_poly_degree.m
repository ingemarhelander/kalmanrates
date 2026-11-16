
% test_select_poly_degree.m
% Testar att select_poly_degree väljer rätt polynomgrad.

clear;
clc;
close all;

%% Lägg till sökväg till funktionerna

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'scripts', 'line_search'));
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

%% Skapa testdata

rng(42);

m = 40;
alpha = sort(rand(m,1));

% Verkligt tredjegradspolynom
% F_true = 2 + 3*alpha - 1.5*alpha.^2 + 0.7*alpha.^3;

% Tillhörande derivata
% S_true = 3 - 3*alpha + 2.1*alpha.^2;

% Verklig femtegradspolynom
F_true = 2 + 3*alpha - 1.5*alpha.^2 ...
         + 0.7*alpha.^3 + 0.4*alpha.^4 - 1.0*alpha.^5;

% Tillhörande derivata
S_true = 3 - 3*alpha + 2.1*alpha.^2 ...
         + 1.6*alpha.^3 - 5*alpha.^4;

noiseF = 0.1;
noiseG = 0.1;

% Lägg till brus
F = F_true + noiseF*randn(m,1);
S = S_true + noiseG*randn(m,1);

%% Kör mult_q_poly_fits

qMax = 5;

fits = mult_q_poly_fits(alpha, F, S, qMax);

% Metod A: tröskelbaserad urval
[qA, fitA, infoA] = select_poly_degree(fits, 'threshold');

% Metod B: poängbaserad urval
[qB, fitB, infoB] = select_poly_degree(fits, 'score');

%% Skriv ut resultat

fprintf('\nRESULTAT FRÅN POLYNOMVAL\n');

if isempty(qA)
    fprintf('Metod A hittade ingen godkänd grad.\n');
else
    fprintf('Metod A valde grad: %d\n', qA);
end

if isempty(qB)
    fprintf('Metod B hittade ingen giltig grad.\n');
else
    fprintf('Metod B valde grad: %d\n', qB);
end

%% Kontrollera resultat

assert(~isempty(qA), 'Metod A valde ingen grad.');
assert(~isempty(qB), 'Metod B valde ingen grad.');

assert(~isempty(qA), 'Metod A valde ingen grad.');
assert(~isempty(qB), 'Metod B valde ingen grad.');

fprintf('\nAlla tester passerade!\n');

fprintf('\nPoäng från Metod B:\n');

for k = 1:length(infoB.scores)
    fprintf('Grad %d: %.4f\n', k, infoB.scores(k));
end