% demo_lr - Enkel visning av gradval med LR-test för olika antal punkter.
%
% Ett fjärdegradspolynom mäts i m = 2, 3, ..., 10 punkter. För varje m
% väljer select_poly_degree_lr graden, och en figur visar mätpunkterna,
% det valda ML-polynomet, den sanna funktionen och residualerna som streck
% (rött = punkten ovanför kurvan, blått = nedanför).
%
% Titeln visar LR för varje steg q -> q+1. Graden höjs om LR > 3.84.
% Högsta grad är m - 2, så rätt grad (4) kan väljas först från m = 6.
% Sista raden visar vad lag-1-autokorrelationen hade valt.
%
% Byt F_true och S_true för att prova en annan funktion.
% Kör från projektroten med:  demo_lr

clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

F_true = @(a) 1 - 4*a.^3 + 5*a.^4;
S_true = @(a) -12*a.^2 + 20*a.^3;
trueName = 'fjärdegrads';

sF = 0.003;   % brus på F
sS = 0.015;   % brus på S

% Färger som syns på både ljus och mörk bakgrund
red    = [1.0 0.3 0.3];
blue   = [0.3 0.6 1.0];
orange = [1.0 0.6 0.0];
green  = [0.2 0.8 0.2];
gray   = [0.6 0.6 0.6];

alphaPlot = linspace(0, 1, 300)';
fprintf('%3s  %6s  %6s\n', 'm', 'LR', 'lag-1');

for m = 2:10
    rng(m);
    alpha = linspace(0, 1, m)';
    F = F_true(alpha) + sF*randn(m, 1);
    S = S_true(alpha) + sS*randn(m, 1);

    q = [];  fit = [];  info = [];  qLag = [];
    evalc('[q, fit, info] = select_poly_degree_lr(alpha, F, S);');   % utan utskrifter
    evalc('qLag = select_poly_degree(alpha, F, S);');
    fprintf('%3d  %6d  %6d\n', m, q, qLag);

    % LR för varje steg q -> q+1
    steps = strings(0);
    for k = 1:numel(info.LR)
        qFrom = 1 + k;
        if info.LR(k) > info.crit
            steps(end+1) = sprintf('%d\\rightarrow%d: LR = %.1f > %.2f, höj', qFrom, qFrom + 1, info.LR(k), info.crit); %#ok<SAGROW>
        else
            steps(end+1) = sprintf('%d\\rightarrow%d: LR = %.1f \\leq %.2f, stanna', qFrom, qFrom + 1, info.LR(k), info.crit); %#ok<SAGROW>
        end
    end
    qMax = max(2, min(6, m - 2));
    if q == qMax
        steps(end+1) = sprintf('högsta grad m - 2 = %d nådd', qMax);
    end

    % Figur
    figure('Name', sprintf('m = %d', m));
    hold on;
    p = polyval(flipud(fit.y), alpha);
    for j = 1:m
        if fit.rF(j) >= 0
            color = red;
        else
            color = blue;
        end
        plot([alpha(j) alpha(j)], [p(j) F(j)], 'Color', color, 'LineWidth', 2);
    end
    hTrue = plot(alphaPlot, F_true(alphaPlot), '--', 'Color', gray, 'LineWidth', 1);
    hFit  = plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), '-', 'Color', green, 'LineWidth', 2);
    hData = plot(alpha, F, 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 7);
    grid on;
    hold off;
    xlabel('Steglängd \alpha');
    ylabel('F');
    legend([hData hFit hTrue], {'Mätningar', sprintf('ML-polynom, q = %d', q), ...
            sprintf('Sann funktion (%s)', trueName)}, 'Location', 'best');
    title({sprintf('m = %d punkter:  LR-testet väljer q = %d', m, q), ...
           char(strjoin(steps, ',   ')), ...
           sprintf('(lag-1-autokorrelation hade valt q = %d)', qLag)});
end
