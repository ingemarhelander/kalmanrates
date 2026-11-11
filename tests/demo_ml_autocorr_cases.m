% demo_ml_autocorr_cases - Flera fall av gradval med ML och autokorrelation.
%
% Samma idé som demo_ml_autocorr, men för flera påhittade funktioner.
% Varje fall ger en figur med en ruta per provad grad (show_ml_degrees).
% I kommandofönstret skrivs vilka grader som provades.
%
% Förväntat (enligt simulering, kan skilja med slumpen):
%   1. Kvadratisk          2             rätt direkt, residualerna är brus
%   2. Kubisk              2 -> 3
%   3. Fjärdegrads, 16 p   2 -> 3 -> 4   fler punkter, tydligare mönster
%   4. sin(3a)             2 -> 3 -> 4   inget polynom, grad 4 räcker
%   5. 1 - exp(-4a)        2 -> 3 -> 4   planar ut, som nära ett maximum
%   6. Fjärdegrads, 8 p    2 -> 3        för få punkter: grad 4 hittas inte
%   7. Ojämna punkter      2 -> 3        punkterna klumpade som i en line search
%
% Kör från projektroten med:  demo_ml_autocorr_cases

clear;
clc;
close all;

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));
addpath(testDir);

sF = 0.003;   % brus på F
sS = 0.015;   % brus på S

% Sanna funktioner F och derivatan S
quad    = {@(a) 1 + 3*a - 5*a.^2,      @(a) 3 - 10*a};
cubic   = {@(a) 1 + 2*a - 6*a.^3,      @(a) 2 - 18*a.^2};
quartic = {@(a) 1 - 4*a.^3 + 5*a.^4,   @(a) -12*a.^2 + 20*a.^3};
sine    = {@(a) sin(3*a),              @(a) 3*cos(3*a)};
expo    = {@(a) 1 - exp(-4*a),         @(a) 4*exp(-4*a)};

% Fall: namn, funktion, steglängder, brus på F och S
cases = {
    '1. Kvadratisk (kontroll), 12 punkter', quad,    linspace(0, 1, 12)', sF,    sS
    '2. Kubisk, 12 punkter',                cubic,   linspace(0, 1, 12)', sF,    sS
    '3. Fjärdegrads, 16 punkter',           quartic, linspace(0, 1, 16)', sF,    sS
    '4. sin(3\alpha), 12 punkter',          sine,    linspace(0, 1, 12)', 0.001, 0.003
    '5. 1 - exp(-4\alpha), 12 punkter',     expo,    linspace(0, 1, 12)', 0.001, 0.004
    '6. Fjärdegrads, bara 8 punkter',       quartic, linspace(0, 1, 8)',  sF,    sS
    '7. Fjärdegrads, ojämna punkter',       quartic, [0 0.05 0.1 0.2 0.4 0.5 0.55 0.6 0.62 0.65 0.8 1]', sF, sS
};

for c = 1:size(cases, 1)
    [name, fun, alpha, noiseF, noiseS] = cases{c, :};
    rng(c);
    F = fun{1}(alpha) + noiseF*randn(size(alpha));
    S = fun{2}(alpha) + noiseS*randn(size(alpha));
    qs = [];
    evalc('qs = show_ml_degrees(alpha, F, S, fun{1}, name);');   % utan utskrifter från poly_fit_ml
    fprintf('%-40s grader: %s\n', name, strjoin(string(qs), ' -> '));
end
