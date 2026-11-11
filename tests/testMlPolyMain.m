classdef testMlPolyMain < matlab.unittest.TestCase
% testMlPolyMain - Tester för ML-boxen (ml_poly_main) med lag-1-gradval:
% att ML-skattningen i poly_fit_ml och lag1_autocorr fungerar tillsammans.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testMlPolyMain.m')

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
            % Olika brus på F och S ska skattas, och q = 2 väljas.
            sF = 0.01;  sS = 0.1;
            alpha = linspace(0, 1, 30)';
            q = zeros(20, 1);  ratioF = q;  ratioS = q;
            for i = 1:20
                rng(i);
                [F, S] = quadratic(alpha, sF, sS);
                [q(i), fit] = mlQuiet(alpha, F, S);
                ratioF(i) = fit.sigmaF/sF;
                ratioS(i) = fit.sigmaG/sS;
            end
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.8);
            testCase.verifyEqual(median(ratioF), 1, 'AbsTol', 0.3);
            testCase.verifyEqual(median(ratioS), 1, 'AbsTol', 0.3);
        end

        function falseAlarmRateIsNearLevel(testCase)
            % Rätt grad och bara brus: graden ska sällan höjas.
            for m = [6 20]
                alpha = linspace(0, 1, m)';
                raised = false(300, 1);
                for i = 1:300
                    rng(i);
                    [F, S] = quadratic(alpha, 0.01, 0.1);
                    raised(i) = mlQuiet(alpha, F, S) > 2;
                end
                testCase.verifyLessThanOrEqual(mean(raised), 0.12, sprintf('m = %d', m));
            end
        end

        function tooLowDegreeIsRaised(testCase)
            % sin är inget polynom: q = 2 räcker inte.
            alpha = linspace(0, 1, 12)';
            q = zeros(50, 1);
            for i = 1:50
                rng(i);
                F = sin(3*alpha) + 1e-3*randn(size(alpha));
                S = 3*cos(3*alpha) + 1e-3*randn(size(alpha));
                q(i) = mlQuiet(alpha, F, S);
            end
            testCase.verifyGreaterThanOrEqual(mean(q > 2), 0.9);
        end

        function fewPointsDoNotCrash(testCase)
            % 2-6 punkter: giltig grad, ändliga koefficienter, konvergens.
            for m = 2:6
                rng(m);
                alpha = sort(rand(m, 1));
                [F, S] = quadratic(alpha, 0.01, 0.1);
                [q, fit] = mlQuiet(alpha, F, S);
                testCase.verifyGreaterThanOrEqual(q, 2, sprintf('m = %d', m));
                testCase.verifyLessThanOrEqual(q, 6, sprintf('m = %d', m));
                testCase.verifyTrue(all(isfinite(fit.y)), sprintf('m = %d', m));
                testCase.verifyTrue(fit.converged, sprintf('m = %d', m));
            end
        end

        function exactPolynomialGivesDegreeTwo(testCase)
            % Utan brus: q = 2, rätt koefficienter, konvergens.
            alpha = linspace(0, 1, 10)';
            [F, S] = quadratic(alpha, 0, 0);
            [q, fit] = mlQuiet(alpha, F, S);
            testCase.verifyEqual(q, 2);
            testCase.verifyTrue(fit.converged);
            testCase.verifyEqual(fit.y, [1; 3; -5], 'AbsTol', 1e-8);
        end

        function noiseOnlyOnS(testCase)
            % F exakt och brus bara på S: ska inte ge för hög grad.
            alpha = linspace(0, 1, 15)';
            q = zeros(20, 1);
            for i = 1:20
                rng(i);
                [F, S] = quadratic(alpha, 0, 0.1);
                [q(i), fit] = mlQuiet(alpha, F, S);
                testCase.verifyTrue(all(isfinite(fit.y)), sprintf('frö %d', i));
            end
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.8);
        end
    end
end

function [F, S] = quadratic(alpha, sF, sS)
% 1 + 3*alpha - 5*alpha^2 och derivatan, med brus.
    F = 1 + 3*alpha - 5*alpha.^2 + sF*randn(size(alpha));
    S = 3 - 10*alpha + sS*randn(size(alpha));
end

function [q, fit] = mlQuiet(alpha, F, S)
% ml_poly_main med lag-1, utan utskrifterna från poly_fit_ml.
    q = [];  fit = [];
    evalc('[q, fit] = ml_poly_main(alpha, F, S);');
end
