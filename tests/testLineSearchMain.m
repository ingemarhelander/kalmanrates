classdef testLineSearchMain < matlab.unittest.TestCase
% testLineSearchMain - Tester för hela line_search_main på en påhittad
% målfunktion. tests/stubs/kalmanObj.m ersätter Kalman.
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
        function findsBetterPointWithBothMethods(testCase)
            for method = {'lag1', 'lr'}
                rng(1);
                K = objective(0.05, 0.1);
                [f0, g0] = K([0; 0]);
                [alpha, xNew, fNew] = runQuiet([0; 0], f0, g0, [2; 0.5], K, method{1});
                testCase.verifyGreaterThan(alpha, 0, method{1});
                testCase.verifyGreaterThan(fNew, f0, method{1});
                testCase.verifyEqual(xNew, [0; 0] + alpha*[2; 0.5], 'AbsTol', 1e-12, method{1});
            end
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

        function minPointsGivesMorePoints(testCase)
            % Med minPoints = 6 får sökningen inte stoppa på lutning eller
            % brus förrän 6 punkter finns, bara om nästa punkt är för nära.
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, 'lr', struct('minPoints', 6));
            testCase.verifyTrue(numel(evals.A) >= 6 || contains(evals.reason, 'för nära'), evals.reason);
        end

        function nonFinitePointsAreSkipped(testCase)
            % f är NaN för långa steg. Sökningen ska inte krascha och ska
            % returnera en ändlig punkt (Req. 42).
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
    end
end

function K = objective(sF, sG)
% f(x) = h(x1) - 2*x2^2 med brus, och gradienten.
    h  = @(t) 6*(1 - exp(-3*t)) - 2.5*t.^2;
    dh = @(t) 18*exp(-3*t) - 5*t;
    K = @(x) deal(h(x(1)) - 2*x(2)^2 + sF*randn, [dh(x(1)); -4*x(2)] + sG*randn(2, 1));
end

function [f, g] = nanForLongSteps(K, x)
% Som K, men NaN när x1 > 1 (motsvarar alpha > 0.5 längs d = [2; 0.5]).
    [f, g] = K(x);
    if x(1) > 1
        f = NaN;
    end
end

function [alpha, xNew, fNew, gNew, evals] = runQuiet(x, f, g, d, K, method, lsOpts)
% line_search_main utan utskrifter.
    if nargin < 7
        lsOpts = [];
    end
    alpha = [];  xNew = [];  fNew = [];  gNew = [];  evals = [];
    evalc('[alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, [], [], method, lsOpts);');
end
