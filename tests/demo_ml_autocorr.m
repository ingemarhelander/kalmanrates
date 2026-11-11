% demo_ml_autocorr - Enkel visning av ML-polynom och lag-1-autokorrelation.
%
% Sann funktion är ett fjärdegradspolynom. Vi anpassar ML-polynom av grad
% q = 2, 3, ... och ritar en figur per grad: mätpunkterna, det anpassade
% polynomet och residualerna som streck från kurvan till punkterna
% (rött = punkten ovanför kurvan, blått = nedanför).
%
% För låg grad: punkterna ligger i klumpar ovanför eller nedanför kurvan
% (+ + + - - - + +), grannar har samma tecken och rho är stort.
% Rätt grad: tecknen hoppar slumpmässigt (+ - + + - + - -) och rho är nära 0.
%
% Kör från projektroten med:  demo_ml_autocorr

clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

% Mockdata: 10 punkter av ett fjärdegradspolynom med brus
rng(1);
m = 10;
alpha = linspace(0, 1, m)';

F_true = @(a) 1 - 2*a + 12*a.^2 - 20*a.^3 + 10*a.^4;
S_true = @(a) -2 + 24*a - 60*a.^2 + 40*a.^3;

F = F_true(alpha) + 0.01*randn(m, 1);
S = S_true(alpha) + 0.05*randn(m, 1);

alphaPlot = linspace(0, 1, 300)';

for q = 2:6
    fit = poly_fit_ml(alpha, F, S, q, false);
    rho = lag1_autocorr([fit.rF, fit.rG]);
    tau = whiteness_threshold(alpha, q, fit);
    ok = max(rho) <= tau;

    p = polyval(flipud(fit.y), alpha);          % polynomet i mätpunkterna
    signs = repmat('-', 1, m);
    signs(fit.rF >= 0) = '+';

    figure('Name', sprintf('q = %d', q));
    hold on;

    % Residualer: streck från kurvan till punkten
    for j = 1:m
        if fit.rF(j) >= 0
            color = 'r';
        else
            color = 'b';
        end
        plot([alpha(j) alpha(j)], [p(j) F(j)], color, 'LineWidth', 2);
    end

    hTrue = plot(alphaPlot, F_true(alphaPlot), 'k--', 'LineWidth', 1);
    hFit  = plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), 'g-', 'LineWidth', 2);
    hData = plot(alpha, F, 'ko', 'MarkerFaceColor', 'k');

    xlabel('Steglängd \alpha');
    ylabel('F');
    legend([hData hFit hTrue], {'Mätningar', sprintf('ML-polynom, q = %d', q), ...
            'Sann funktion (grad 4)'}, 'Location', 'best');
    grid on;
    hold off;

    if ok
        verdict = 'OK, residualerna ser ut som brus';
    else
        verdict = 'höj graden';
    end
    title({sprintf('q = %d:   tecken på residualerna:   %s', q, strjoin(cellstr(signs'), ' ')), ...
           sprintf('\\rho_F = %.2f,  \\rho_S = %.2f,  \\tau = %.2f   =>   %s', rho(1), rho(2), tau, verdict)});

    if ok
        break
    end
end
