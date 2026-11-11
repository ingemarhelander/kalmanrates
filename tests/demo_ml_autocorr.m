% demo_ml_autocorr - Enkel visning av ML-polynom och lag-1-autokorrelation.
%
% Sann funktion är ett fjärdegradspolynom. Vi anpassar ML-polynom av grad
% q = 2, 3, ... och ritar en figur per grad: mätpunkterna, det anpassade
% polynomet och residualerna som streck från kurvan till punkterna
% (rött = punkten ovanför kurvan, blått = nedanför).
% Om graden höjs för långt (5 eller 6) är det slumpen: byt frö i rng.
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

% Mockdata: 12 punkter av ett fjärdegradspolynom med brus.
% Med färre punkter syns mönstret för dåligt: residualerna för grad q byter
% tecken minst q+1 gånger, så med t.ex. 10 punkter blir klumparna så korta
% att rho inte skiljer sig från brus och graden aldrig höjs.
rng(1);
m = 12;
alpha = linspace(0, 1, m)';

F_true = @(a) 1 - 4*a.^3 + 5*a.^4;
S_true = @(a) -12*a.^2 + 20*a.^3;

F = F_true(alpha) + 0.003*randn(m, 1);
S = S_true(alpha) + 0.015*randn(m, 1);

% Färger som syns på både ljus och mörk bakgrund
red    = [1.0 0.3 0.3];
blue   = [0.3 0.6 1.0];
orange = [1.0 0.6 0.0];
green  = [0.2 0.8 0.2];
gray   = [0.6 0.6 0.6];

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
            color = red;
        else
            color = blue;
        end
        plot([alpha(j) alpha(j)], [p(j) F(j)], 'Color', color, 'LineWidth', 2);
    end

    hTrue = plot(alphaPlot, F_true(alphaPlot), '--', 'Color', gray, 'LineWidth', 1);
    hFit  = plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), '-', 'Color', green, 'LineWidth', 2);
    hData = plot(alpha, F, 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 7);

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
