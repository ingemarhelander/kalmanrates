% demo_few_points - Gradval med få punkter (högst 10): autokorrelation mot LR-test.
%
% Del 1: lag-1-autokorrelation (select_poly_degree, nuvarande metod).
% Del 2: likelihoodkvottest (select_poly_degree_lr, prototyp).
% Båda körs på samma data. En figur per metod: rader = funktioner,
% kolumner = antal punkter m. Varje ruta visar vald grad med mätpunkter,
% polynom, sann funktion och residualer (rött ovanför, blått nedanför).
% Del 3: tabell över 50 slumpade dataset per fall.
%
% Förväntat (enligt simulering):
%   Autokorrelation: väljer nästan alltid q = 2, ser inte att graden är för låg.
%   LR-test: hittar rätt grad från ca 5-6 punkter, falsklarm ca 5 %.
%   Med m punkter är högsta grad i LR-testet m - 2 (t.ex. m = 5 => q <= 3).
%
% Kör från projektroten med:  demo_few_points

clear;
clc;
close all;

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));

sF = 0.003;   % brus på F
sS = 0.015;   % brus på S
ms = [4 5 6 8 10];

% Namn, sann F, sann S, sann grad (NaN = inget polynom)
funs = {
    'Kvadratisk',    @(a) 1 + 3*a - 5*a.^2,    @(a) 3 - 10*a,             2
    'Kubisk',        @(a) 1 + 2*a - 6*a.^3,    @(a) 2 - 18*a.^2,          3
    'Fjärdegrads',   @(a) 1 - 4*a.^3 + 5*a.^4, @(a) -12*a.^2 + 20*a.^3,   4
    '1 - exp(-4a)',  @(a) 1 - exp(-4*a),       @(a) 4*exp(-4*a),          NaN
};
selectors = {'Autokorrelation (lag-1)', @select_poly_degree
           'LR-test (prototyp)',      @select_poly_degree_lr};

%% Del 1 och 2: en figur per metod
for k = 1:size(selectors, 1)
    figure('Name', selectors{k, 1}, 'Position', [30 30 1600 900]);
    t = tiledlayout(size(funs, 1), numel(ms), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('%s: vald grad för få punkter', selectors{k, 1}));
    for i = 1:size(funs, 1)
        for j = 1:numel(ms)
            m = ms(j);
            alpha = linspace(0, 1, m)';
            rng(100*i + m);
            F = funs{i, 2}(alpha) + sF*randn(m, 1);
            S = funs{i, 3}(alpha) + sS*randn(m, 1);
            [q, fit, info] = quietSelect(selectors{k, 2}, alpha, F, S);
            nexttile(t);
            plotTile(alpha, F, funs{i, 2}, fit);
            title(sprintf('%s, m = %d:  q = %s', funs{i, 1}, m, strjoin(string(info.q), '\rightarrow')));
        end
    end
end

%% Del 3: tabell över många dataset
nTrials = 50;
fprintf('Andel dataset (av %d) där rätt grad valdes. För 1 - exp(-4a): medel av vald grad.\n\n', nTrials);
fprintf('%-14s %4s   %-10s %-10s\n', 'Funktion', 'm', 'Autokorr', 'LR-test');
for i = 1:size(funs, 1)
    for j = 1:numel(ms)
        m = ms(j);
        alpha = linspace(0, 1, m)';
        qq = zeros(nTrials, 2);
        for s = 1:nTrials
            rng(s);
            F = funs{i, 2}(alpha) + sF*randn(m, 1);
            S = funs{i, 3}(alpha) + sS*randn(m, 1);
            for k = 1:2
                qq(s, k) = quietSelect(selectors{k, 2}, alpha, F, S);
            end
        end
        qTrue = funs{i, 4};
        if isnan(qTrue)
            fprintf('%-14s %4d   q = %-6.1f q = %-6.1f\n', funs{i, 1}, m, mean(qq(:, 1)), mean(qq(:, 2)));
        else
            fprintf('%-14s %4d   %5.0f %%    %5.0f %%\n', funs{i, 1}, m, 100*mean(qq == qTrue));
        end
    end
    fprintf('\n');
end
fprintf('Kvadratisk: andelen som INTE är rätt är falsklarm (borde vara ca 5 %%).\n');
fprintf('LR-test: högsta grad är m - 2, så fjärdegrads kan inte hittas med m = 4 eller 5.\n');

%% Lokala funktioner

function [q, fit, info] = quietSelect(fun, alpha, F, S)
% Kör en gradvalsfunktion utan utskrifterna från poly_fit_ml.
    q = [];  fit = [];  info = [];
    evalc('[q, fit, info] = fun(alpha, F, S);');
end

function plotTile(alpha, F, F_true, fit)
% Mätpunkter, polynom, sann funktion och residualer. Färger syns på både
% ljus och mörk bakgrund.
    red = [1.0 0.3 0.3];  blue = [0.3 0.6 1.0];
    orange = [1.0 0.6 0.0];  green = [0.2 0.8 0.2];  gray = [0.6 0.6 0.6];
    p = polyval(flipud(fit.y), alpha);
    hold on;
    for j = 1:numel(alpha)
        if fit.rF(j) >= 0, c = red; else, c = blue; end
        plot([alpha(j) alpha(j)], [p(j) F(j)], 'Color', c, 'LineWidth', 2);
    end
    a = linspace(0, 1, 200)';
    plot(a, F_true(a), '--', 'Color', gray);
    plot(a, polyval(flipud(fit.y), a), '-', 'Color', green, 'LineWidth', 1.5);
    plot(alpha, F, 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 5);
    grid on;
    hold off;
end
