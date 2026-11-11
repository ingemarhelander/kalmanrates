classdef testMlPolyAnimations < matlab.unittest.TestCase
% testMlPolyAnimations - Animeringar som visar hur ML-boxen arbetar.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testMlPolyAnimations.m')
% Eller ett test i taget, t.ex.
%   runtests('tests/testMlPolyAnimations.m', 'ProcedureName', 'degreeByDegree')
%
%   degreeByDegree   grad för grad: polynomet, residualerna och lag-plotten
%                    när graden höjs, tills residualerna ser ut som brus.
%                    Avslutas med en översikt med alla grader bredvid varandra.
%   noiseVsStructure histogram som växer fram: rho för rent brus mot rho
%                    när graden är för låg, med tröskeln tau.
%   noiseSweep       samma funktion med mer och mer brus: bruset döljer
%                    strukturen och vald grad sjunker.
%   lineSearchMovie  hela line_search_main på en påhittad målfunktion
%                    (attrapp för kalmanObj), varv för varv.
%
% Varje test kontrollerar också något, så det går att köra med
% ShowPlots = false som vanliga tester.

    properties (Constant)
        ShowPlots = true;   % rita figurer och animeringar
        PauseTime = 1.0;    % sekunder per bild
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'scripts', 'line_search')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'tests')));           % plotMlBox
        end
    end

    methods (TestMethodSetup)
        function restoreRandomState(testCase)
            state = rng;
            testCase.addTeardown(@() rng(state));
        end
    end

    methods (Test)
        function degreeByDegree(testCase)
            % sin(3*alpha) med litet brus kräver hög grad. Anpassa varje grad
            % q = 2..6 och visa hur rho sjunker under tau när graden räcker.
            alpha = linspace(0, 1, 20)';
            sTrue = 1e-3;
            rng(2);
            F = sin(3*alpha) + sTrue*randn(size(alpha));
            S = 3*cos(3*alpha) + sTrue*randn(size(alpha));

            qs = 2:6;  n = numel(qs);
            fits = cell(n, 1);  rho = zeros(n, 2);  tau = zeros(n, 1);
            fprintf('\n  degreeByDegree (sin(3*alpha), m = 20, brus %.0e):\n', sTrue);
            fprintf('     q   rho_F   rho_S     tau   RMS r_F   beslut\n');
            for k = 1:n
                fits{k} = fitQuiet(alpha, F, S, qs(k));
                rho(k, :) = lag1_autocorr([fits{k}.rF, fits{k}.rG]);
                tau(k) = whiteness_threshold(alpha, qs(k), fits{k});
                fprintf('    %2d  %6.2f  %6.2f  %6.2f  %8.1e   %s\n', qs(k), rho(k, 1), rho(k, 2), ...
                        tau(k), rmsOf(fits{k}.rF), verdict(max(rho(k, :)) <= tau(k)));
            end
            qChosen = mlQuiet(alpha, F, S);
            fprintf('    ml_poly_main valde q = %d\n', qChosen);

            if testMlPolyAnimations.ShowPlots
                animateDegrees(alpha, F, S, qs, fits, rho, tau, sTrue, qChosen);
                filmStrip(alpha, F, S, qs, fits, rho, tau, sTrue, qChosen);
            end

            % q = 2 ska förkastas, och ml_poly_main ska välja den första grad
            % där rho är under tau (samma beslut som i tabellen).
            kOk = find(max(rho, [], 2) <= tau, 1);
            testCase.verifyGreaterThan(max(rho(1, :)), tau(1));
            testCase.verifyNotEmpty(kOk);
            testCase.verifyEqual(qChosen, qs(kOk));
        end

        function noiseVsStructure(testCase)
            % Histogram som växer fram. Blått: rätt grad (q = 2 på en
            % kvadratisk funktion), bara brus kvar. Rött: för låg grad
            % (q = 2 på sin), struktur kvar. tau ska skilja dem åt.
            alpha = linspace(0, 1, 10)';
            nTrials = 400;
            rhoW = nan(nTrials, 1);  rhoS = rhoW;  tauW = rhoW;  tauS = rhoW;
            if testMlPolyAnimations.ShowPlots
                figure('Name', 'Brus mot struktur', 'Position', [100 100 900 500]);
                edges = -1:0.05:1;
                hW = histogram(rhoW, edges, 'FaceColor', [0.2 0.4 0.85], 'FaceAlpha', 0.6);
                hold on;  grid on;
                hS = histogram(rhoS, edges, 'FaceColor', [0.85 0.2 0.2], 'FaceAlpha', 0.6);
                hT = xline(0, 'k-', '\tau', 'LineWidth', 2);
                xlabel('max(\rho_F, \rho_S) vid q = 2');  ylabel('antal');
                legend([hW hS], {'rätt grad: bara brus', 'för låg grad: sin(3\alpha)'}, 'Location', 'northwest');
            end
            for i = 1:nTrials
                rng(i);
                Fw = 1 + 3*alpha - 5*alpha.^2 + 0.01*randn(size(alpha));
                Sw = 3 - 10*alpha + 0.1*randn(size(alpha));
                fit = fitQuiet(alpha, Fw, Sw, 2);
                rhoW(i) = max(lag1_autocorr([fit.rF, fit.rG]));
                tauW(i) = whiteness_threshold(alpha, 2, fit);

                Fs = sin(3*alpha) + 0.01*randn(size(alpha));
                Ss = 3*cos(3*alpha) + 0.01*randn(size(alpha));
                fit = fitQuiet(alpha, Fs, Ss, 2);
                rhoS(i) = max(lag1_autocorr([fit.rF, fit.rG]));
                tauS(i) = whiteness_threshold(alpha, 2, fit);

                if testMlPolyAnimations.ShowPlots && (mod(i, 10) == 0 || i == nTrials)
                    hW.Data = rhoW(1:i);  hS.Data = rhoS(1:i);
                    hT.Value = median([tauW(1:i); tauS(1:i)]);
                    title(sprintf(['%d försök. Över \\tau: rätt grad %.1f %% (falsklarm), ' ...
                                   'för låg grad %.1f %% (upptäckt)'], i, ...
                                  100*mean(rhoW(1:i) > tauW(1:i)), 100*mean(rhoS(1:i) > tauS(1:i))));
                    drawnow;
                end
            end
            falseAlarm = mean(rhoW > tauW);
            power = mean(rhoS > tauS);
            fprintf('\n  noiseVsStructure: falsklarm %.1f %%, upptäckt %.1f %%\n', 100*falseAlarm, 100*power);
            testCase.verifyLessThanOrEqual(falseAlarm, 0.12);
            testCase.verifyGreaterThan(power, falseAlarm + 0.2);
        end

        function noiseSweep(testCase)
            % Samma funktion och samma slumptal, men bruset skalas upp. Med
            % lite brus syns att sin inte är ett andragradspolynom, med mycket
            % brus går det inte att se och q = 2 väljs.
            alpha = linspace(0, 1, 15)';
            rng(1);
            eF = randn(size(alpha));  eS = randn(size(alpha));
            sigmas = logspace(-4, -0.5, 12);
            q = zeros(size(sigmas));
            fig = [];
            for k = 1:numel(sigmas)
                F = sin(3*alpha) + sigmas(k)*eF;
                S = 3*cos(3*alpha) + sigmas(k)*eS;
                q(k) = mlQuiet(alpha, F, S);
                if testMlPolyAnimations.ShowPlots
                    fig = plotMlBox(alpha, F, S, sprintf('Brusnivå %.1e', sigmas(k)), fig);
                    pause(testMlPolyAnimations.PauseTime);
                end
            end
            fprintf('\n  noiseSweep: vald grad per brusnivå: %s\n', mat2str(q));
            if testMlPolyAnimations.ShowPlots
                figure('Name', 'Vald grad mot brus');
                semilogx(sigmas, q, 'o-', 'LineWidth', 2, 'MarkerFaceColor', 'auto');
                grid on;  yticks(2:6);  ylim([1.5 6.5]);
                xlabel('brusnivå \sigma');  ylabel('vald grad q');
                title('Mer brus döljer strukturen: graden sjunker');
            end
            testCase.verifyGreaterThan(q(1), 2);
            testCase.verifyLessThanOrEqual(q(end), q(1));
        end

        function lineSearchMovie(testCase)
            % Kör den riktiga line_search_main på en påhittad målfunktion i
            % två dimensioner. tests/stubs/kalmanObj.m ersätter Kalman.
            % Efteråt spelas varven upp: punkterna, ML-polynomet, vald grad
            % och nästa steglängd.
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'tests', 'stubs')));

            c = [0.5; 0.3];        % maximum utan bumpen
            sF = 0.05;  sG = 0.05; % brus på f och g
            fTrue = @(x) -50*sum((x - c).^2, 1) + 2*cos(6*x(1, :));
            gTrue = @(x) -100*(x - c) + [-12*sin(6*x(1, :)); zeros(1, size(x, 2))];
            K = @(x) deal(fTrue(x) + sF*randn, gTrue(x) + sG*randn(2, 1));

            rng(1);
            x0 = [0; 0];
            [f0, g0] = K(x0);
            d = c/0.6;             % sökriktning, maximum ungefär vid alpha = 0.6
            fprintf('\n  lineSearchMovie:\n');
            [alpha, ~, fNew, ~, evals] = line_search_main(x0, f0, g0, d, K, [], []);

            alphas = (d'*(evals.X - x0))/(d'*d);   % steglängd för varje punkt
            S = d'*evals.G;
            n = numel(alphas);
            fprintf('    %d punkter, returnerad alpha = %.3f, f - f0 = %+.2f\n', n, alpha, fNew - f0);

            if testMlPolyAnimations.ShowPlots
                replayLineSearch(alphas, evals.F, S, alpha, @(a) fTrue(x0 + d*a), ...
                                 @(a) d'*gTrue(x0 + d*a));
            end
            testCase.verifyTrue(isfinite(fNew));
            testCase.verifyGreaterThanOrEqual(alpha, 0);
            testCase.verifyLessThanOrEqual(n, 21);
            testCase.verifyGreaterThan(fNew, f0);
        end
    end
end

%% ----------------------------- Animeringar ------------------------------

function animateDegrees(alpha, F, S, qs, fits, rho, tau, sTrue, qChosen)
% En bild per grad. Övre raden F, nedre raden S. Kolumner: anpassning,
% residualer (röd = positiv, blå = negativ), lag-plot r_{j+1} mot r_j.
% Höger: rho och tau för graderna hittills.
    fig = figure('Name', 'Grad för grad', 'Position', [60 80 1400 680]);
    aa = linspace(min(alpha), max(alpha), 300)';
    n = numel(qs);
    for k = 1:n
        clf(fig);
        t = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        ok = max(rho(k, :)) <= tau(k);
        title(t, sprintf('q = %d:   max(\\rho_F, \\rho_S) = %.2f  %s  \\tau = %.2f   =>   %s', ...
              qs(k), max(rho(k, :)), ternary(ok, '\leq', '>'), tau(k), verdict(ok)), ...
              'Color', ternary(ok, [0 0.5 0], [0.75 0 0]), 'FontWeight', 'bold');

        % F
        ax = nexttile(t, 1);  hold(ax, 'on');  grid(ax, 'on');
        for j = 1:k - 1
            plot(ax, aa, poly_value(fits{j}, aa), 'Color', [0.8 0.8 0.8]);
        end
        plot(ax, aa, poly_value(fits{k}, aa), 'r', 'LineWidth', 2);
        plot(ax, alpha, F, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        title(ax, 'F och ML-polynom (grått = lägre grader)');  xlabel(ax, '\alpha');
        residualPlot(nexttile(t, 2), alpha, fits{k}.rF, sTrue, 'r_F');
        lagPlot(nexttile(t, 3), fits{k}.rF, rho(k, 1), tau(k), '\rho_F');

        % S
        ax = nexttile(t, 5);  hold(ax, 'on');  grid(ax, 'on');
        for j = 1:k - 1
            plot(ax, aa, polyDeriv(fits{j}, aa), 'Color', [0.8 0.8 0.8]);
        end
        plot(ax, aa, polyDeriv(fits{k}, aa), 'r', 'LineWidth', 2);
        plot(ax, alpha, S, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
        title(ax, 'S = g''d och derivatan');  xlabel(ax, '\alpha');
        residualPlot(nexttile(t, 6), alpha, fits{k}.rG, sTrue, 'r_S');
        lagPlot(nexttile(t, 7), fits{k}.rG, rho(k, 2), tau(k), '\rho_S');

        % rho per grad
        ax = nexttile(t, 4, [2 1]);  hold(ax, 'on');  grid(ax, 'on');
        plot(ax, qs(1:k), tau(1:k), 'k--d', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
        plot(ax, qs(1:k), rho(1:k, 1), '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
        plot(ax, qs(1:k), rho(1:k, 2), '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
        yline(ax, 0, 'k:');
        if qs(k) >= qChosen
            xline(ax, qChosen, 'g-', 'vald', 'LineWidth', 2);
        end
        xlim(ax, [qs(1) - 0.5, qs(end) + 0.5]);  ylim(ax, [-1 1]);  xticks(ax, qs);
        xlabel(ax, 'grad q');
        title(ax, 'Autokorrelation per grad');
        legend(ax, {'\tau', '\rho_F', '\rho_S'}, 'Location', 'southwest');

        drawnow;
        pause(2*testMlPolyAnimations.PauseTime);
    end
end

function filmStrip(alpha, F, S, qs, fits, rho, tau, sTrue, qChosen)
% Alla grader bredvid varandra: anpassning, residualer för F och S, lag-plot.
    n = numel(qs);
    figure('Name', 'Översikt: alla grader', 'Position', [40 40 300*n 900]);
    t = tiledlayout(4, n, 'TileSpacing', 'tight', 'Padding', 'compact');
    title(t, 'Översikt: när graden höjs försvinner mönstret i residualerna');
    aa = linspace(min(alpha), max(alpha), 300)';
    for k = 1:n
        ok = max(rho(k, :)) <= tau(k);
        ax = nexttile(t, k);  hold(ax, 'on');  grid(ax, 'on');
        plot(ax, aa, poly_value(fits{k}, aa), 'r', 'LineWidth', 1.5);
        plot(ax, alpha, F, 'k.', 'MarkerSize', 10);
        label = sprintf('q = %d   (\\rho = %.2f, \\tau = %.2f)\n%s', qs(k), max(rho(k, :)), tau(k), verdict(ok));
        if qs(k) == qChosen
            label = [label, '  <- VALD'];
        end
        title(ax, label, 'Color', ternary(ok, [0 0.5 0], [0.75 0 0]));
        residualPlot(nexttile(t, n + k), alpha, fits{k}.rF, sTrue, 'r_F');
        residualPlot(nexttile(t, 2*n + k), alpha, fits{k}.rG, sTrue, 'r_S');
        ax = nexttile(t, 3*n + k);
        lagPlot(ax, fits{k}.rF, rho(k, 1), tau(k), '\rho_F');
        lagPlot(ax, fits{k}.rG, rho(k, 2), tau(k), '\rho_S', [0.85 0.45 0.1]);
        title(ax, sprintf('\\rho_F = %.2f, \\rho_S = %.2f', rho(k, 1), rho(k, 2)));
    end
end

function replayLineSearch(alphas, F, S, alphaRet, fLine, sLine)
% Spelar upp line search varv för varv. Varv k: ML-boxen har sett de k
% första punkterna och nästa punkt är k + 1.
    fig = figure('Name', 'Line search varv för varv', 'Position', [80 80 1300 650]);
    n = numel(alphas);
    aa = linspace(0, max(alphas)*1.1, 300);
    ok = isfinite(F(:)') & isfinite(S(:)');
    for k = 2:n
        clf(fig);
        t = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        idx = find(ok(1:k));
        last = k == n;
        q = 2;  fit = [];  info = [];
        if numel(idx) >= 2
            [q, fit, info] = mlQuiet(alphas(idx), F(idx), S(idx));
        end
        if last
            heading = sprintf('Klar efter %d punkter: returnerar \\alpha = %.3f', n, alphaRet);
        else
            heading = sprintf('Varv %d: %d punkter, vald grad q = %d, nästa \\alpha = %.3f', ...
                              k - 1, k, q, alphas(k + 1));
        end
        title(t, heading, 'FontWeight', 'bold');

        % f längs linjen
        ax = nexttile(t, 1, [1 2]);  hold(ax, 'on');  grid(ax, 'on');
        plot(ax, aa, fLine(aa), '--', 'Color', [0.6 0.6 0.6], 'DisplayName', 'sann f (utan brus)');
        if ~isempty(fit)
            plot(ax, aa, poly_value(fit, aa), 'r', 'LineWidth', 2, 'DisplayName', sprintf('ML-polynom, q = %d', q));
        end
        plot(ax, alphas(1:k), F(1:k), 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'utvärderade punkter');
        text(ax, alphas(1:k), F(1:k), arrayfun(@(i) sprintf('  %d', i), 0:k - 1, 'UniformOutput', false), 'FontSize', 9);
        if ~last
            xline(ax, alphas(k + 1), 'b-', 'nästa \alpha', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        else
            [~, j] = min(abs(alphas - alphaRet));
            plot(ax, alphas(j), F(j), 'gp', 'MarkerSize', 22, 'MarkerFaceColor', 'g', 'DisplayName', 'returnerad punkt');
        end
        xlabel(ax, '\alpha');  ylabel(ax, 'f(x + \alpha d)');
        title(ax, 'Målfunktionen längs sökriktningen (siffra = ordning)');
        legend(ax, 'Location', 'southeast');

        % rho per provad grad
        ax = nexttile(t, 3, [2 1]);  hold(ax, 'on');  grid(ax, 'on');
        if ~isempty(info)
            plot(ax, info.q, info.tau, 'k--d', 'MarkerFaceColor', 'k', 'LineWidth', 1.5);
            plot(ax, info.q, info.rho(:, 1), '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
            plot(ax, info.q, info.rho(:, 2), '-s', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
            legend(ax, {'\tau', '\rho_F', '\rho_S'}, 'Location', 'southwest');
        end
        yline(ax, 0, 'k:');
        xlim(ax, [1.5 6.5]);  ylim(ax, [-1 1]);  xticks(ax, 2:6);
        xlabel(ax, 'provad grad q');
        title(ax, ternary(sum(ok(1:k)) < 4, 'Färre än 4 punkter: bara q = 2', 'Autokorrelation per provad grad'));

        % S längs linjen
        ax = nexttile(t, 4, [1 2]);  hold(ax, 'on');  grid(ax, 'on');
        plot(ax, aa, sLine(aa), '--', 'Color', [0.6 0.6 0.6]);
        if ~isempty(fit)
            plot(ax, aa, polyDeriv(fit, aa(:)), 'r', 'LineWidth', 2);
        end
        plot(ax, alphas(1:k), S(1:k), 'ko', 'MarkerFaceColor', 'k');
        yline(ax, 0, 'k:');
        xlabel(ax, '\alpha');  ylabel(ax, 'S = g''d');
        title(ax, 'Lutning längs sökriktningen');

        drawnow;
        pause(1.5*testMlPolyAnimations.PauseTime);
    end
end

%% ----------------------------- Hjälpfunktioner --------------------------

function residualPlot(ax, x, r, sTrue, name)
% Residualer som staplar, röd = positiv och blå = negativ. Grönt band = ±2
% gånger det sanna bruset. Långa rader med samma färg = autokorrelation.
    hold(ax, 'on');  grid(ax, 'on');
    fill(ax, [x(1) x(end) x(end) x(1)], 2*sTrue*[-1 -1 1 1], [0.85 0.95 0.85], 'EdgeColor', 'none');
    pos = r >= 0;
    stem(ax, x(pos), r(pos), 'filled', 'Color', [0.85 0.2 0.2], 'MarkerSize', 4);
    stem(ax, x(~pos), r(~pos), 'filled', 'Color', [0.2 0.4 0.85], 'MarkerSize', 4);
    yline(ax, 0, 'k-');
    lim = max(1.2*max(abs(r)), 4*sTrue);
    ylim(ax, [-lim lim]);  xlim(ax, [x(1) x(end)]);
    title(ax, sprintf('%s (RMS %.1e, grönt = \\pm2\\sigma brus)', name, rmsOf(r)));
end

function lagPlot(ax, r, rho, tau, name, color)
% r_{j+1} mot r_j, normerat. Punkter längs diagonalen = grannar har samma
% tecken (autokorrelation). Ett runt moln = brus. Linjen har lutning rho.
    if nargin < 6, color = [0.2 0.4 0.85]; end
    hold(ax, 'on');  grid(ax, 'on');
    z = r/max(std(r), realmin);
    scatter(ax, z(1:end-1), z(2:end), 25, color, 'filled', 'MarkerFaceAlpha', 0.7);
    plot(ax, [-3 3], rho*[-3 3], '-', 'Color', color, 'LineWidth', 2);
    xlim(ax, [-3 3]);  ylim(ax, [-3 3]);  axis(ax, 'square');
    xlabel(ax, 'r_j');  ylabel(ax, 'r_{j+1}');
    title(ax, sprintf('Lag-plot: %s = %.2f (\\tau = %.2f)', name, rho, tau));
end

function s = polyDeriv(fit, a)
% Derivatan av polynomet. Antar att fit.y är koefficienter för 1, alpha, ...
    y = fit.y;
    dy = (1:numel(y) - 1)'.*y(2:end);
    s = a(:).^(0:numel(dy) - 1)*dy;
end

function v = rmsOf(r)
    v = sqrt(mean(r(:).^2));
end

function s = verdict(ok)
    s = ternary(ok, 'OK: ser ut som brus', 'höj graden: mönster kvar');
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end

function fit = fitQuiet(alpha, F, S, q)
% poly_fit_ml för en given grad, utan figurer och utskrifter.
    fit = [];
    evalc('fit = poly_fit_ml(alpha, F, S, q, false);');
end

function [q, fit, info] = mlQuiet(alpha, F, S)
% ml_poly_main utan utskrifterna från poly_fit_ml.
    q = [];  fit = [];  info = [];
    evalc('[q, fit, info] = ml_poly_main(alpha, F, S);');
end
