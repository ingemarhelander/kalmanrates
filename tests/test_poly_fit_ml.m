
% test_poly_fit_ml - Ritar ML-polynom av olika grad för slumpad testdata.
% Kör från projektroten med:  test_poly_fit_ml
clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

% Generate random mock data
m = 40;

% Randomly distributed step lengths
alpha = sort(0.02 + 0.48*rand(m,1));

% True underlying polynomial
F_true = 1 + 2*alpha - 1.5*alpha.^2 + 0.8*alpha.^3;

% True derivative
S_true = 2 - 3*alpha + 2.4*alpha.^2;

% Add random measurement noise
sigmaF = 0.06;
sigmaG = 0.15;

F = F_true + sigmaF*randn(m,1);
S = S_true + sigmaG*randn(m,1);

% Test polynomial degrees 1 to 5
qMax = 10;
fit = cell(qMax,1);

for q = 1:qMax
    fit{q} = poly_fit_ml(alpha, F, S, q, false);
end

% Points for smooth curves
alphaPlot = linspace(min(alpha), max(alpha), 300)';

% Plot all polynomial fits
figure;

% Function values
subplot(2,1,1);
plot(alpha, F, 'ko', 'MarkerFaceColor', 'k');
hold on;

for q = 1:qMax
    Fplot = polyval(flipud(fit{q}.y), alphaPlot);
    plot(alphaPlot, Fplot, 'LineWidth', 1.5);
end

xlabel('Step length \alpha');
ylabel('Function value F');
title('ML polynomial fits - function values');

legend(['Measurements', arrayfun(@(q) ...
    sprintf('q = %d',q), 1:qMax, 'UniformOutput', false)], ...
    'Location', 'best');

grid on;
hold off;

% Derivative values
subplot(2,1,2);
plot(alpha, S, 'ko', 'MarkerFaceColor', 'k');
hold on;

for q = 1:qMax
    y = fit{q}.y;
    dy = (1:q)' .* y(2:end);

    Splot = polyval(flipud(dy), alphaPlot);
    plot(alphaPlot, Splot, 'LineWidth', 1.5);
end

xlabel('Step length \alpha');
ylabel('Derivative S');
title('ML polynomial fits - derivatives');

legend(['Measurements', arrayfun(@(q) ...
    sprintf('q = %d',q), 1:qMax, 'UniformOutput', false)], ...
    'Location', 'best');

grid on;
hold off;
