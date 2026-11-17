classdef testLineSearchMain < matlab.unittest.TestCase
% testLineSearchMain - Tester för hela line_search_main på en påhittad
% målfunktion. tests/stubs/kalmanObj.m ersätter Kalman.
%
% Målfunktionen längs d = [2; 0.5] har sitt maximum vid alpha ~ 0.29.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testLineSearchMain.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(testDir), 'scripts', 'line_search')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testDir, 'stubs')));
        end
    end

    methods (TestMethodSetup)
        function restoreRandomState(testCase)
            state = rng;
            testCase.addTeardown(@() rng(state));
        end
    end

    methods (Test)
        function findsBetterPointWithAllMethods(testCase)
            for method = {'lag1', 'lr', 'score', 'threshold'}
                rng(1);
                K = objective(0.05, 0.1);
                [f0, g0] = K([0; 0]);
                [alpha, xNew, fNew] = runQuiet([0; 0], f0, g0, [2; 0.5], K, method{1});
                testCase.verifyGreaterThan(alpha, 0, method{1});
                testCase.verifyGreaterThan(fNew, f0, method{1});
                testCase.verifyEqual(xNew, [0; 0] + alpha*[2; 0.5], 'AbsTol', 1e-12, method{1});
            end
        end

        function firstPointIsLongestStep(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr');
            testCase.verifyEqual(evals.A(2), 1);
        end

        function evaluatesOnlyPolynomialMaxima(testCase)
            % Efter första punkten ska varje utvärderad alpha vara förslaget
            % från polynomet i varvet innan.
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr');
            n = numel(evals.A);
            testCase.verifyGreaterThan(n, 2);
            testCase.verifyEqual(evals.A(3:n), evals.AHat(2:n-1), 'AbsTol', 1e-12);
        end

        function findsMaximumWithLittleNoise(testCase)
            rng(1);
            K = objective(1e-3, 1e-3);
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr');
            testCase.verifyLessThanOrEqual(abs(alpha - trueMaximum()), 0.02);
            testCase.verifyTrue(contains(evals.reason, 'ligger nära'), evals.reason);
        end

        function neverRepeatsAnAlpha(testCase)
            % Deterministisk målfunktion, som riktiga Kalman: samma alpha
            % ger samma värden. Ingen alpha får utvärderas två gånger.
            for method = {'lr', 'score', 'threshold'}
                K = deterministicObjective();
                [f0, g0] = K([0; 0]);
                [~, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, method{1});
                testCase.verifyEqual(numel(unique(evals.A)), numel(evals.A), method{1});
                testCase.verifyGreaterThan(fNew, f0, method{1});
            end
        end

        function steepStartShrinksUntilImprovement(testCase)
            % Som i riktiga körningen: en varians nära noll ger en barriärterm
            % med startlutning ~1e15, och alpha = 1 är mycket sämre än
            % startpunkten. Maximum ligger vid alpha = 1e-3. Sökningen ska
            % krympa steget med faktor 10 tills F ökar, inte ge upp.
            K = @(x) deal(0.01*log(x(1) + 1e-17) - 10*x(1) - x(2)^2, ...
                          [0.01/(x(1) + 1e-17) - 10; -2*x(2)]);
            [f0, g0] = K([0; 0]);
            [alpha, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, 'lr');
            testCase.verifyEqual(evals.A(2:4), [1 0.1 0.01], 'AbsTol', 1e-15);
            testCase.verifyGreaterThan(alpha, 0);
            testCase.verifyGreaterThan(fNew, f0);
        end

        function neverStopsWithoutImprovementBeforeMaxEvals(testCase)
            % Alla stopp utom maxEvals kräver en punkt bättre än startpunkten.
            for seed = 1:5
                rng(seed);
                K = objective(0.05, 0.1);
                [f0, g0] = K([0; 0]);
                [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr');
                testCase.verifyTrue(alpha > 0 || contains(evals.reason, 'maxEvals'), evals.reason);
            end
        end

        function stopsAtAlphaOneWhenAlwaysIncreasing(testCase)
            % f ökar hela vägen: polynomet har inget maximum före 1, så
            % förslaget blir 1 igen och sökningen stoppar direkt.
            rng(1);
            K = @(x) deal(x(1) + 1e-3*randn, [1; 0] + 1e-3*randn(2, 1));
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, 'lr');
            testCase.verifyEqual(alpha, 1);
            testCase.verifyEqual(numel(evals.A), 2);
        end

        function maxEvalsLimitsEvaluations(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr', ...
                                           struct('maxEvals', 2, 'dMin', 1e-9));
            testCase.verifyEqual(numel(evals.A), 3);
            testCase.verifyTrue(contains(evals.reason, 'maxEvals'), evals.reason);
        end

        function evalsAreConsistent(testCase)
            % evals har en kolumn per punkt och samma längd på A, Q och AHat.
            rng(2);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lag1');
            n = size(evals.X, 2);
            testCase.verifyEqual(numel(evals.F), n);
            testCase.verifyEqual(size(evals.G, 2), n);
            testCase.verifyEqual(numel(evals.A), n);
            testCase.verifyEqual(numel(evals.Q), n);
            testCase.verifyEqual(numel(evals.AHat), n);
            testCase.verifyEqual(evals.A(1), 0);
            testCase.verifyLessThanOrEqual(n, 21);
            testCase.verifyNotEmpty(evals.reason);
            q = evals.Q(~isnan(evals.Q));
            testCase.verifyTrue(all(q >= 2 & q <= 6));
        end

        function nonFinitePointsAreSkipped(testCase)
            % f är NaN för alpha > 0.5, så första punkten (alpha = 1) är NaN.
            % Sökningen ska halvera, inte krascha, och returnera en ändlig
            % punkt (Req. 42).
            rng(3);
            Kok = objective(0.05, 0.1);
            K = @(x) nanForLongSteps(Kok, x);
            [f0, g0] = K([0; 0]);
            [alpha, ~, fNew, gNew, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lag1');
            testCase.verifyTrue(isfinite(fNew));
            testCase.verifyTrue(all(isfinite(gNew)));
            testCase.verifyLessThanOrEqual(alpha, 0.5);
            testCase.verifyTrue(any(~isfinite(evals.F)));
        end

        function oldOptionsGiveError(testCase)
            % Parametrar från den gamla linjesökningen ska inte ignoreras tyst.
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr', struct('minPoints', 5)), ?MException);
        end
    end
end

function K = objective(sF, sG)
% f(x) = h(x1) - 2*x2^2 med brus, och gradienten.
    h  = @(t) 6*(1 - exp(-3*t)) - 2.5*t.^2;
    dh = @(t) 18*exp(-3*t) - 5*t;
    K = @(x) deal(h(x(1)) - 2*x(2)^2 + sF*randn, [dh(x(1)); -4*x(2)] + sG*randn(2, 1));
end

function K = deterministicObjective()
% Som objective, men "bruset" är en fast funktion av x i stället för randn.
    h  = @(t) 6*(1 - exp(-3*t)) - 2.5*t.^2;
    dh = @(t) 18*exp(-3*t) - 5*t;
    K = @(x) deal(h(x(1)) - 2*x(2)^2 + 0.05*sin(997*x(1)), ...
                  [dh(x(1)); -4*x(2)] + 0.1*[cos(991*x(1)); sin(983*x(2))]);
end

function a = trueMaximum()
% Maximum av f längs d = [2; 0.5] utan brus.
    fLine = @(a) 6*(1 - exp(-6*a)) - 2.5*(2*a).^2 - 2*(0.5*a).^2;
    a = fminbnd(@(a) -fLine(a), 0, 1);
end

function [f, g] = nanForLongSteps(K, x)
% Som K, men NaN när x1 > 1 (motsvarar alpha > 0.5 längs d = [2; 0.5]).
    [f, g] = K(x);
    if x(1) > 1
        f = NaN;
    end
end

function [alpha, xNew, fNew, gNew, evals] = runQuiet(x, f, g, d, K, method, opts)
% line_search_main utan utskrifter. method blir opts.degreeMethod.
    if nargin < 7
        opts = struct();
    end
    opts.degreeMethod = method;
    alpha = [];  xNew = [];  fNew = [];  gNew = [];  evals = [];
    evalc('[alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, [], [], opts);');
end
