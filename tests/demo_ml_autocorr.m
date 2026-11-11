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

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

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
