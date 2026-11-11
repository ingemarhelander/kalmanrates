% demo_lr - Gradval med LR-test för m = 2, 3, ..., 10 punkter, en figur per m.
%
% Ett fjärdegradspolynom mäts med brus i m punkter. select_poly_degree_lr
% väljer graden och figuren visar polynomet med residualerna som streck.
% Titeln visar LR för varje steg q -> q+1 och vilken grad lag-1 hade valt.
%
% Kör från projektroten med:  demo_lr

clear;
clc;
close all;

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));
addpath(testDir);

F_true = @(a) 1 - 4*a.^3 + 5*a.^4;
S_true = @(a) -12*a.^2 + 20*a.^3;

fprintf('%3s  %6s  %6s\n', 'm', 'LR', 'lag-1');
for m = 2:10
    rng(m);
    alpha = linspace(0, 1, m)';
    F = F_true(alpha) + 0.003*randn(m, 1);
    S = S_true(alpha) + 0.015*randn(m, 1);

    q = [];  fit = [];  info = [];  qLag = [];
    evalc('[q, fit, info] = ml_poly_main(alpha, F, S, [], ''lr'');');   % utan utskrifter
    evalc('qLag = ml_poly_main(alpha, F, S, [], ''lag1'');');
    fprintf('%3d  %6d  %6d\n', m, q, qLag);

    % LR för varje steg q -> q+1
    steps = {};
    for k = 1:numel(info.LR)
        if info.LR(k) > info.crit
            steps{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f > %.2f, höj', k + 1, k + 2, info.LR(k), info.crit); %#ok<SAGROW>
        else
            steps{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f \\leq %.2f, stanna', k + 1, k + 2, info.LR(k), info.crit); %#ok<SAGROW>
        end
    end
    if q == info.qMax
        steps{end+1} = sprintf('högsta grad m - 2 = %d nådd', info.qMax); %#ok<SAGROW>
    end

    plot_fit_residuals(alpha, F, fit, F_true, ...
        {sprintf('LR-test, m = %d punkter:  väljer q = %d', m, q), ...
         strjoin(steps, ',   '), ...
         sprintf('(lag-1 hade valt q = %d)', qLag)});
end
