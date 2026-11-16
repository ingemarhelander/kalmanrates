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

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

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

%% Lokala funktioner

function fig = plot_fit_residuals(alpha, F, fit, F_true, titleLines)
% plot_fit_residuals - Ritar en ML-anpassning med residualer i en ny figur.
%
% Mätpunkter (orange), ML-polynomet (grönt), sann funktion (grå, streckad)
% och residualerna som streck från kurvan till punkterna: rött = punkten
% ligger ovanför kurvan, blått = nedanför. Färgerna syns på både ljus och
% mörk bakgrund.
%
%   alpha, F    mätningar
%   fit         anpassning från poly_fit_ml (fit.y, fit.rF)
%   F_true      sann funktion (funktionshandtag)
%   titleLines  rubrik, text eller cell med en rad per element

    red    = [1.0 0.3 0.3];
    blue   = [0.3 0.6 1.0];
    orange = [1.0 0.6 0.0];
    green  = [0.2 0.8 0.2];
    gray   = [0.6 0.6 0.6];

    alpha = alpha(:);  F = F(:);
    p = polyval(flipud(fit.y), alpha);       % polynomet i mätpunkterna
    alphaPlot = linspace(min(alpha), max(alpha), 300)';

    fig = figure;
    hold on;
    for j = 1:numel(alpha)
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
    legend([hData hFit hTrue], {'Mätningar', sprintf('ML-polynom, q = %d', numel(fit.y) - 1), ...
            'Sann funktion'}, 'Location', 'best');
    title(titleLines);
end
