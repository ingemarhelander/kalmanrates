% demo_ml_autocorr - Enkel animering av ML-polynom och lag-1-autokorrelation.
%
% Anpassar ML-polynom av grad q = 2, 3, ... till få brusiga punkter.
% För låg grad: residualerna bildar ett mönster (grannar har samma tecken),
% rho är stort och graden höjs. Rätt grad: residualerna ser ut som brus.
% Kör från projektroten med:  demo_ml_autocorr

clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

% Mockdata: 8 punkter av sin(3*alpha) med litet brus
rng(1);
m = 8;
alpha = linspace(0, 1, m)';
F = sin(3*alpha) + 1e-3*randn(m, 1);
S = 3*cos(3*alpha) + 1e-3*randn(m, 1);

alphaPlot = linspace(0, 1, 300)';
figure('Name', 'ML och autokorrelation');

for q = 2:6
    fit = poly_fit_ml(alpha, F, S, q, false);
    rho = lag1_autocorr([fit.rF, fit.rG]);
    tau = whiteness_threshold(alpha, q, fit);
    ok = max(rho) <= tau;

    % Funktionsvärden och ML-polynomet
    subplot(3, 1, 1);
    plot(alpha, F, 'ko', 'MarkerFaceColor', 'k');
    hold on;
    plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), 'r-', 'LineWidth', 1.5);
    hold off;
    grid on;
    ylabel('F');
    legend('Mätningar', 'ML-polynom', 'Location', 'best');
    if ok
        title(sprintf('q = %d: max \\rho = %.2f \\leq \\tau = %.2f  =>  OK, residualerna ser ut som brus', q, max(rho), tau));
    else
        title(sprintf('q = %d: max \\rho = %.2f > \\tau = %.2f  =>  höj graden', q, max(rho), tau));
    end

    % Residualer för F
    subplot(3, 1, 2);
    plot(alpha, fit.rF, 'bo-', 'MarkerFaceColor', 'b');
    yline(0, 'k-');
    grid on;
    ylabel('r_F');
    title(sprintf('Residualer för F, \\rho_F = %.2f', rho(1)));

    % Residualer för S
    subplot(3, 1, 3);
    plot(alpha, fit.rG, 'mo-', 'MarkerFaceColor', 'm');
    yline(0, 'k-');
    grid on;
    xlabel('Steglängd \alpha');
    ylabel('r_S');
    title(sprintf('Residualer för S, \\rho_S = %.2f', rho(2)));

    drawnow;
    pause(2);

    if ok
        break
    end
end
