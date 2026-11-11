% demo_ml_autocorr - Gradval med lag-1-autokorrelation, en figur per grad.
%
% Ett fjärdegradspolynom mäts med brus i 12 punkter. ML-polynom av grad
% q = 2, 3, ... anpassas tills residualerna ser ut som brus (max rho <= tau).
% Varje figur visar polynomet, residualerna som streck och deras tecken.
%
% Kör från projektroten med:  demo_ml_autocorr

clear;
clc;
close all;

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));
addpath(testDir);

F_true = @(a) 1 - 4*a.^3 + 5*a.^4;
S_true = @(a) -12*a.^2 + 20*a.^3;
m = 12;

rng(1);
alpha = linspace(0, 1, m)';
F = F_true(alpha) + 0.003*randn(m, 1);
S = S_true(alpha) + 0.015*randn(m, 1);

for q = 2:6
    fit = [];
    evalc('fit = poly_fit_ml(alpha, F, S, q, false);');   % utan utskrifter
    rho = lag1_autocorr([fit.rF, fit.rG]);
    tau = whiteness_threshold(alpha, q, fit);
    ok = max(rho) <= tau;

    signs = repmat('-', 1, m);
    signs(fit.rF >= 0) = '+';
    if ok
        verdict = sprintf('max \\rho = %.2f \\leq \\tau = %.2f  =>  OK, residualerna ser ut som brus', max(rho), tau);
    else
        verdict = sprintf('max \\rho = %.2f > \\tau = %.2f  =>  höj graden', max(rho), tau);
    end
    plot_fit_residuals(alpha, F, fit, F_true, ...
        {sprintf('Lag-1, m = %d, q = %d:   tecken:  %s', m, q, strjoin(cellstr(signs'), ' ')), ...
         sprintf('\\rho_F = %.2f,  \\rho_S = %.2f,   %s', rho(1), rho(2), verdict)});

    if ok
        break
    end
end
