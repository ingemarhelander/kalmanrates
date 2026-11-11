classdef testMlPolyMain < matlab.unittest.TestCase
% testMlPolyMain - Testar ML-boxen (ml_poly_main): att gradvalet med
% lag1_autocorr fungerar tillsammans med ML-skattningen i poly_fit_ml.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testMlPolyMain.m')
% Kräver varken data eller mex. Datan är påhittade funktioner med brus.
% Varje test skriver ut en kort sammanfattning.

    methods (TestClassSetup)
        function addPaths(testCase)
            projectRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(projectRoot, 'scripts', 'line_search')));
        end
    end

    methods (TestMethodSetup)
        function restoreRandomState(testCase)
            state = rng;
            testCase.addTeardown(@() rng(state));
        end
    end

    methods (Test)
        function differentNoiseIsEstimated(testCase)
            % Kvadratisk funktion, 10 gånger mer brus på S än på F.
            % ML ska skatta båda brusnivåerna och oftast välja q = 2.
            sF = 0.01;  sS = 0.1;
            alpha = linspace(0, 1, 30)';
            seeds = 1:20;
            q = zeros(size(seeds));  ratioF = q;  ratioS = q;
            for i = seeds
                rng(i);
                [F, S] = quadratic(alpha, sF, sS);
                [q(i), fit] = mlQuiet(alpha, F, S);
                ratioF(i) = fit.sigmaF/sF;
                ratioS(i) = fit.sigmaG/sS;
            end
            fprintf('\n  differentNoiseIsEstimated: q = 2 i %.0f %%, sigmaF/sant = %.2f, sigmaS/sant = %.2f (median)\n', ...
                    100*mean(q == 2), median(ratioF), median(ratioS));
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.8);
            testCase.verifyEqual(median(ratioF), 1, 'AbsTol', 0.3);
            testCase.verifyEqual(median(ratioS), 1, 'AbsTol', 0.3);
        end

        function falseAlarmRateIsNearLevel(testCase)
            % Rätt grad och bara brus: graden ska höjas i ungefär 5 % av
            % fallen (level = 0.95 i select_poly_degree). Testar att tau, som
            % räknas med ML-vikterna, är rätt kalibrerad.
            % m = 6 liknar en riktig linjesökning, m = 20 är fler punkter.
            for m = [6 20]
                alpha = linspace(0, 1, m)';
                nTrials = 300;
                raised = false(nTrials, 1);
                for i = 1:nTrials
                    rng(i);
                    [F, S] = quadratic(alpha, 0.01, 0.1);
                    raised(i) = mlQuiet(alpha, F, S) > 2;
                end
                fprintf('\n  falseAlarmRateIsNearLevel: m = %2d, graden höjd i %.1f %% (förväntat ca 5 %%)\n', ...
                        m, 100*mean(raised));
                testCase.verifyLessThanOrEqual(mean(raised), 0.12, sprintf('m = %d', m));
            end
        end

        function tooLowDegreeIsRaised(testCase)
            % sin är inte ett polynom: q = 2 lämnar ett slätt mönster i
            % residualerna, så graden ska höjas.
            alpha = linspace(0, 1, 12)';
            seeds = 1:50;
            q = zeros(size(seeds));  rho2 = q;
            for i = seeds
                rng(i);
                F = sin(3*alpha) + 1e-3*randn(size(alpha));
                S = 3*cos(3*alpha) + 1e-3*randn(size(alpha));
                q(i) = mlQuiet(alpha, F, S);
                fit2 = quietFit(alpha, F, S, 2);
                rho2(i) = max(lag1_autocorr([fit2.rF, fit2.rG]));
            end
            fprintf('\n  tooLowDegreeIsRaised: q > 2 i %.0f %%, grader %s, rho vid q = 2: %.2f (median)\n', ...
                    100*mean(q > 2), mat2str(unique(q)), median(rho2));
            testCase.verifyGreaterThanOrEqual(mean(q > 2), 0.9);
        end

        function fewPointsDoNotCrash(testCase)
            % Få punkter som i en riktig linjesökning. Ska ge en giltig
            % anpassning, och ML-iterationen ska konvergera.
            fprintf('\n  fewPointsDoNotCrash:\n');
            for m = 2:6
                rng(m);
                alpha = sort(rand(m, 1));
                [F, S] = quadratic(alpha, 0.01, 0.1);
                [q, fit] = mlQuiet(alpha, F, S);
                fprintf('    m = %d: q = %d, iterationer = %4d, sigmaF = %.2e, sigmaS = %.2e\n', ...
                        m, q, fit.iterations, fit.sigmaF, fit.sigmaG);
                testCase.verifyGreaterThanOrEqual(q, 2, sprintf('m = %d', m));
                testCase.verifyLessThanOrEqual(q, 6, sprintf('m = %d', m));
                testCase.verifyTrue(all(isfinite(fit.y)), sprintf('m = %d', m));
                testCase.verifyTrue(fit.converged, sprintf('m = %d', m));
            end
        end

        function exactPolynomialGivesDegreeTwo(testCase)
            % Utan brus passar q = 2 exakt. Ska stanna direkt och konvergera.
            alpha = linspace(0, 1, 10)';
            [F, S] = quadratic(alpha, 0, 0);
            [q, fit] = mlQuiet(alpha, F, S);
            fprintf('\n  exactPolynomialGivesDegreeTwo: q = %d, iterationer = %d\n', q, fit.iterations);
            testCase.verifyEqual(q, 2);
            testCase.verifyTrue(fit.converged);
            testCase.verifyEqual(fit.y, [1; 3; -5], 'AbsTol', 1e-8);
        end

        function noiseOnlyOnS(testCase)
            % F exakt, brus bara på S. sigmaF hamnar vid golvet, vilket inte
            % får göra att anpassningen eller gradvalet går sönder.
            alpha = linspace(0, 1, 15)';
            seeds = 1:20;
            q = zeros(size(seeds));
            for i = seeds
                rng(i);
                [F, S] = quadratic(alpha, 0, 0.1);
                [q(i), fit] = mlQuiet(alpha, F, S);
                testCase.verifyTrue(all(isfinite(fit.y)), sprintf('frö %d', i));
            end
            fprintf('\n  noiseOnlyOnS: q = 2 i %.0f %%, grader %s\n', 100*mean(q == 2), mat2str(unique(q)));
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.8);
        end
    end
end

function [F, S] = quadratic(alpha, sF, sS)
% Kvadratisk funktion 1 + 3*alpha - 5*alpha^2 och dess derivata, med brus.
    F = 1 + 3*alpha - 5*alpha.^2 + sF*randn(size(alpha));
    S = 3 - 10*alpha + sS*randn(size(alpha));
end

function [q, fit] = mlQuiet(alpha, F, S, qPrev)
% ml_poly_main utan utskrifterna från poly_fit_ml.
    if nargin < 4
        qPrev = [];
    end
    q = [];  fit = [];
    evalc('[q, fit] = ml_poly_main(alpha, F, S, qPrev);');
end

function fit = quietFit(alpha, F, S, q)
% poly_fit_ml för en given grad, utan figurer och utskrifter.
    fit = [];
    evalc('fit = poly_fit_ml(alpha, F, S, q, false);');
end
