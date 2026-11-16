
% test_mult_q_poly_fits.m
% Testar att mult_q_poly_fits genererar och sparar
% ML-anpassningar och residualdiagnostik för alla grader.

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

noiseF = 0.01;
noiseG = 0.01;

% Lägg till brus
F = F_true + noiseF*randn(m,1);
S = S_true + noiseG*randn(m,1);

%% Kör mult_q_poly_fits

qMax = 5;

fits = mult_q_poly_fits(alpha, F, S, qMax);

%% Kontrollera antal anpassningar

assert(numel(fits) == qMax, ...
    'Incorrect number of polynomial fits.');

%% Kontrollera att diagnostiken finns

requiredFields = {
    'y'
    'rF'
    'rG'
    'sigmaF'
    'sigmaG'
    'rhoF'
    'rhoG'
    'trendCorrF'
    'trendCorrG'
    'pTrendF'
    'pTrendG'
    'iterations'
    'converged'
};

for q = 1:qMax

    fprintf('\nChecking polynomial degree %d...\n', q);

    for k = 1:numel(requiredFields)
        assert(isfield(fits{q}, requiredFields{k}), ...
            'Missing field: %s in degree %d.', ...
            requiredFields{k}, q);
    end

    assert(fits{q}.converged, ...
        'ML did not converge for degree %d.', q);

    assert(isfinite(fits{q}.rhoF));
    assert(isfinite(fits{q}.rhoG));

    assert(isfinite(fits{q}.pTrendF));
    assert(isfinite(fits{q}.pTrendG));

    fprintf('Degree %d passed!\n', q);

end

fprintf('\nAll mult_q_poly_fits tests passed!\n');
