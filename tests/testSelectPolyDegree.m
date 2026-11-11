classdef testSelectPolyDegree < matlab.unittest.TestCase
% testSelectPolyDegree - Tester för lag1_autocorr och select_poly_degree.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testSelectPolyDegree.m')
% Kräver varken data eller mex. poly_fit_ml skattar vikterna wF och wG med ML.

    methods (TestClassSetup)
        function addPaths(testCase)
            % setupPaths lägger inte till scripts/line_search.
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
        function lag1KnownValues(testCase)
            % Alternerande följd: -3/4. Linjär följd: centrerad [-1.5 -0.5 0.5 1.5]
            % ger (0.75 - 0.25 + 0.75)/5 = 1/4. Kolumnerna räknas var för sig.
            r = [1 1; -1 2; 1 3; -1 4];
            testCase.verifyEqual(lag1_autocorr(r), [-0.75, 0.25], 'AbsTol', 1e-14);
            testCase.verifyEqual(lag1_autocorr(r(:, 1)'), -0.75, 'AbsTol', 1e-14);
        end

        function lag1ConstantSeriesIsWhite(testCase)
            testCase.verifyEqual(lag1_autocorr(ones(5, 1)), 0);
        end

        function quadraticGivesDegreeTwo(testCase)
            alpha = linspace(0, 1, 8)';
            q = select_poly_degree(alpha, 1 - 5*(alpha - 0.3).^2, -10*(alpha - 0.3));
            testCase.verifyEqual(q, 2);
        end

        function quarticRaisesDegree(testCase)
            % Andragradsanpassningen ger släta residualer, alltså stort rho.
            alpha = linspace(0, 1, 15)';
            q = select_poly_degree(alpha, alpha - 2*alpha.^4, 1 - 8*alpha.^3);
            testCase.verifyGreaterThan(q, 2);
        end

        function fewPointsGiveDegreeTwo(testCase)
            % Färre än fyra punkter: rho beräknas inte (Req. 28).
            for m = 2:3
                alpha = linspace(0, 1, m)';
                q = select_poly_degree(alpha, alpha - 2*alpha.^4, 1 - 8*alpha.^3);
                testCase.verifyEqual(q, 2, sprintf('m = %d', m));
            end
        end

        function noiseAroundQuadraticMostlyKeepsDegreeTwo(testCase)
            % Med rätt grad ska testet sällan slå larm. Nivån 95 % ger i
            % storleksordningen 5 % falsklarm, så kravet 80 % har marginal.
            alpha = linspace(0, 1, 10)';
            seeds = 1:50;
            q = zeros(size(seeds));
            for i = seeds
                rng(i);
                F = 1 - 5*(alpha - 0.3).^2 + 0.05*randn(size(alpha));
                S = -10*(alpha - 0.3) + 0.05*randn(size(alpha));
                q(i) = select_poly_degree(alpha, F, S);
            end
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.8);
        end

        function orderOfPointsDoesNotMatter(testCase)
            alpha = linspace(0, 1, 15)';
            F = alpha - 2*alpha.^4;
            S = 1 - 8*alpha.^3;
            [q1, fit1] = select_poly_degree(alpha, F, S);
            p = [15 3 9 1 12 6 2 14 8 5 11 4 13 7 10];
            [q2, fit2] = select_poly_degree(alpha(p), F(p), S(p));
            testCase.verifyEqual(q2, q1);
            testCase.verifyEqual(fit2.y, fit1.y, 'AbsTol', 1e-12);
        end

        function hysteresisNeverLowersDegree(testCase)
            % Med förra varvets grad qPrev = 4 är det svårare att välja en
            % lägre grad. Graden kan alltså bara bli lika eller högre (Req. 30).
            alpha = linspace(0, 1, 10)';
            for seed = 1:20
                rng(seed);
                F = 1 - 5*(alpha - 0.3).^2 + 0.05*randn(size(alpha));
                S = -10*(alpha - 0.3) + 0.05*randn(size(alpha));
                qWithout = select_poly_degree(alpha, F, S, []);
                qWith = select_poly_degree(alpha, F, S, 4);
                testCase.verifyGreaterThanOrEqual(qWith, qWithout, sprintf('frö %d', seed));
            end
        end

        function globalRandomStateIsUnchanged(testCase)
            % Gradvalet får inte påverka slumptalen, annars ger samma frö
            % olika resultat i linjesökningen.
            rng(3);
            before = rng;
            alpha = linspace(0, 1, 10)';
            select_poly_degree(alpha, alpha - 2*alpha.^4, 1 - 8*alpha.^3);
            testCase.verifyEqual(rng, before);
        end
    end
end
