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
        function findsBetterPoint(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [alpha, xNew, fNew] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            testCase.verifyGreaterThan(alpha, 0);
            testCase.verifyGreaterThan(fNew, f0);
            testCase.verifyEqual(xNew, [0; 0] + alpha*[2; 0.5], 'AbsTol', 1e-12);
        end

        function firstPointIsLongestStep(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            testCase.verifyEqual(evals.A(2), 1);
        end

        function evaluatesOnlyPolynomialMaxima(testCase)
            % Efter första punkten ska varje utvärderad alpha vara förslaget
            % från polynomet i varvet innan.
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            n = numel(evals.A);
            testCase.verifyGreaterThan(n, 2);
            testCase.verifyEqual(evals.A(3:n), evals.AHat(2:n-1), 'AbsTol', 1e-12);
        end

        function findsMaximumWithLittleNoise(testCase)
            rng(1);
            K = objective(1e-3, 1e-3);
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            testCase.verifyLessThanOrEqual(abs(alpha - trueMaximum()), 0.02);
            % Stopp på liten lutning, i bästa punkten eller vid polynomets max
            testCase.verifyTrue(contains(evals.reason, 'ligger nära') || ...
                                contains(evals.reason, 'är liten'), evals.reason);
        end

        function neverRepeatsAnAlpha(testCase)
            % Deterministisk målfunktion, som riktiga Kalman: samma alpha
            % ger samma värden. Ingen alpha får utvärderas två gånger.
            K = deterministicObjective();
            [f0, g0] = K([0; 0]);
            [~, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            testCase.verifyEqual(numel(unique(evals.A)), numel(evals.A));
            testCase.verifyGreaterThan(fNew, f0);
        end

        function steepStartShrinksUntilImprovement(testCase)
            % Som i riktiga körningen: en varians nära noll ger en barriärterm
            % med startlutning ~1e15, och alpha = 1 är mycket sämre än
            % startpunkten. Maximum ligger vid alpha = 1e-3. Sökningen ska
            % halvera steget tills F ökar, inte ge upp.
            K = @(x) deal(0.01*log(x(1) + 1e-17) - 10*x(1) - x(2)^2, ...
                          [0.01/(x(1) + 1e-17) - 10; -2*x(2)]);
            [f0, g0] = K([0; 0]);
            [alpha, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K);
            testCase.verifyEqual(evals.A(2:4), [1 0.5 0.25], 'AbsTol', 1e-15);
            testCase.verifyGreaterThan(alpha, 0);
            testCase.verifyGreaterThan(fNew, f0);
        end

        function neverStopsWithoutImprovementBeforeMaxEvals(testCase)
            % Alla stopp utom maxEvals kräver en punkt bättre än startpunkten.
            for seed = 1:5
                rng(seed);
                K = objective(0.05, 0.1);
                [f0, g0] = K([0; 0]);
                [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
                testCase.verifyTrue(alpha > 0 || contains(evals.reason, 'maxEvals'), evals.reason);
            end
        end

        function stopsAtAlphaOneWhenAlwaysIncreasing(testCase)
            % f ökar hela vägen: polynomet har inget maximum före 1, så
            % förslaget blir 1 igen och sökningen stoppar direkt.
            rng(1);
            K = @(x) deal(x(1) + 1e-3*randn, [1; 0] + 1e-3*randn(2, 1));
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K);
            testCase.verifyEqual(alpha, 1);
            testCase.verifyEqual(numel(evals.A), 2);
        end

        function maxEvalsLimitsEvaluations(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, ...
                                           struct('maxEvals', 2, 'dMin', 1e-9));
            testCase.verifyEqual(numel(evals.A), 3);
            testCase.verifyTrue(contains(evals.reason, 'maxEvals'), evals.reason);
        end

        function evalsAreConsistent(testCase)
            % evals har en kolumn per punkt och samma längd på A, Q och AHat.
            rng(2);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
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
            [alpha, ~, fNew, gNew, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K);
            testCase.verifyTrue(isfinite(fNew));
            testCase.verifyTrue(all(isfinite(gNew)));
            testCase.verifyLessThanOrEqual(alpha, 0.5);
            % NaN-punkten (alpha = 1) är utvärderad men inte med i evals
            testCase.verifyTrue(all(isfinite(evals.F)) && all(isfinite(evals.G(:))));
            testCase.verifyFalse(any(evals.A == 1));
        end

        function firstStepFromInitialStepLength(testCase)
            % Första steget är deltaF/(g'*d) = 100/S0, begränsat till
            % [alphaMin, 1].
            K = @(x) deal(-0.01*(x(1) - 30)^2, [-0.02*(x(1) - 30); 0]);
            [f0, g0] = K([0; 0]);              % S0 = 0.6 längs d = [1; 0]
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                                           struct('alphaMin', 0.1));
            testCase.verifyEqual(evals.A(2), 1);    % 100/0.6 > 1
            K = @(x) deal(-1e3*(x(1) - 1)^2, [-2e3*(x(1) - 1); 0]);
            [f0, g0] = K([0; 0]);              % S0 = 2000
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                                           struct('alphaMin', 1e-3));
            testCase.verifyEqual(evals.A(2), 0.05, 'AbsTol', 1e-15);  % 100/2000
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                                           struct('alphaMin', 0.1));
            testCase.verifyEqual(evals.A(2), 0.1);  % 100/2000 < alphaMin
        end

        function complexValuesAreInvalid(testCase)
            % Ett för långt steg ger log av ett negativt tal (gränsen ligger
            % vid alpha = 0.5). Sådana punkter ska sparas som NaN, inte bli
            % förbättringar. Startlutning 3, så första steget blir 1.
            Kd = @(x) deal(log(0.5 - x(1)) - 0.1*x(1)^2 + 5*x(1), [-1/(0.5 - x(1)) - 0.2*x(1) + 5; 0]);
            [f0, g0] = Kd([0; 0]);
            [~, ~, fNew, gNew, evals] = runQuiet([0; 0], f0, g0, [1; 0], Kd);
            testCase.verifyTrue(isreal(fNew) && isreal(gNew));
            testCase.verifyGreaterThan(fNew, f0);
            % alpha = 1 gav komplext F och är borttagen ur evals
            testCase.verifyFalse(any(evals.A == 1));
            testCase.verifyTrue(isreal(evals.F) && all(isfinite(evals.F)));
        end

        function nonAscentDirectionEvaluatesNothing(testCase)
            % dF/dalpha(0) < 0: ingen Kalman-körning, startpunkten returneras.
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [alpha, xNew, fNew, ~, evals] = runQuiet([0; 0], f0, g0, -[2; 0.5], K);
            testCase.verifyEqual(alpha, 0);
            testCase.verifyEqual(xNew, [0; 0]);
            testCase.verifyEqual(fNew, f0);
            testCase.verifyEqual(numel(evals.A), 1);
            testCase.verifyTrue(contains(evals.reason, 'uppåtriktning'), evals.reason);
        end

        function repeatedInvalidProposalGoesHalfway(testCase)
            % F ökar fram till maximum vid alpha ~ 0.8, och gränsen ligger
            % vid alpha = 1 (log(1 - alpha) = -Inf där). Efter krympningen
            % till 0.1 föreslår polynomet gränsen igen. Sökningen ska gå
            % halvvägs mot den i stället för att stoppa vid 0.1.
            Kd = @(x) deal(log(1 - x(1)) - 0.1*x(1)^2 + 5*x(1), [-1/(1 - x(1)) - 0.2*x(1) + 5; 0]);
            [f0, g0] = Kd([0; 0]);
            [alpha, ~, fNew] = runQuiet([0; 0], f0, g0, [1; 0], Kd);
            testCase.verifyGreaterThan(alpha, 0.5);
            testCase.verifyGreaterThan(fNew, f0);
        end

        function neverStepsPastAlphaMax(testCase)
            % F ökar fram till alpha = 30, men steget får aldrig bli längre
            % än alphaMax (1 i runQuiet, eller 5).
            K = @(x) deal(-0.01*(x(1) - 30)^2, [-0.02*(x(1) - 30); 0]);
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                                               struct('alphaMin', 0.1));
            testCase.verifyLessThanOrEqual(max(evals.A), 1);
            testCase.verifyEqual(alpha, 1);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                                               struct('alphaMin', 0.1, 'alphaMax', 5));
            testCase.verifyEqual(evals.A(2), 5);    % 100/0.6 > 5
            testCase.verifyLessThanOrEqual(max(evals.A), 5);
            testCase.verifyEqual(alpha, 5);
        end

        function halvingFollowsSlopeAtBest(testCase)
            % Varje halvvägs-steg ska gå åt det håll som lutningen i den
            % bästa punkten pekar: positiv lutning -> större alpha.
            Kb = @(x) deal(10*x(1) + 0.01*log(1 - x(1)), [10 - 0.01/(1 - x(1)); 0]);
            [f0, g0] = Kb([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], Kb, ...
                                           struct('alphaMin', 0.1, 'dMin', 1e-3));
            % Lutningen är positiv ända fram till 0.999, så alla giltiga
            % punkter efter den första ska ligga till höger om den bästa före
            A = evals.A(2:end);
            S = evals.G(1, 2:end);
            for j = 2:numel(A)
                [~, b] = max(evals.F(2:j));        % bästa punkten före punkt j
                if S(b) > 0 && evals.F(j+1) < max(evals.F(2:j)) && A(j) < A(b)
                    testCase.verifyFail(sprintf('punkt %d (alpha = %.4g) ligger till vänster om bästa', j, A(j)));
                end
            end
        end

        function halvingStopsAtBoundary(testCase)
            % F ökar fram till maximum vid alpha = 0.999, strax före gränsen
            % alpha = 1 (log(1 - alpha) ger -Inf/komplext F där). Sökningen
            % ska stoppa nära maximum, antingen för att intervallet är mindre
            % än dMin (inringat) eller för att lutningen där är liten, och
            % inte gå till maxEvals.
            Kb = @(x) deal(10*x(1) + 0.01*log(1 - x(1)), [10 - 0.01/(1 - x(1)); 0]);
            [f0, g0] = Kb([0; 0]);
            [alpha, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], Kb, ...
                                                  struct('alphaMin', 0.1, 'dMin', 1e-3));
            testCase.verifyTrue(contains(evals.reason, 'inringat') || ...
                                contains(evals.reason, 'är liten'), evals.reason);
            testCase.verifyGreaterThan(alpha, 0.99);
            testCase.verifyGreaterThan(fNew, f0);
        end

        function oldOptionsGiveError(testCase)
            % Parametrar från den gamla linjesökningen ska inte ignoreras tyst.
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, struct('minPoints', 5)), ?MException);
            % xLbd finns inte längre: d är redan skalad mot variansgränserna.
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, struct('xLbd', [0; 0])), ?MException);
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

function [alpha, xNew, fNew, gNew, evals] = runQuiet(x, f, g, d, K, opts)
% line_search_main utan utskrifter.
% Om inget annat anges är alphaMin = alphaMax = 1, så första steget blir
% alpha = 1 och testerna ovan inte beror på initial_step_length.
    if nargin < 6
        opts = struct();
    end
    if ~isfield(opts, 'alphaMax'), opts.alphaMax = 1; end
    if ~isfield(opts, 'alphaMin'), opts.alphaMin = 1; end
    alpha = [];  xNew = [];  fNew = [];  gNew = [];  evals = [];
    evalc('[alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, [], [], opts);');
end
