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
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(testDir), 'scripts', 'functions')));
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
            % 0.05: de förvalda punkterna ligger jämnt i [0, alpha1], så
            % sökningen kan stoppa när polynomets max hamnar nära en av dem
            % i stället för precis vid toppen 0.294.
            testCase.verifyLessThanOrEqual(abs(alpha - trueMaximum()), 0.05);
            % Stopp på liten lutning, i bästa punkten eller vid polynomets max
            testCase.verifyTrue(contains(evals.reason, 'ligger nära') || ...
                contains(evals.reason, 'är liten'), evals.reason);
        end

        function newtonFindsMaximumWithLittleNoise(testCase)
            % Samma som ovan med ls_maxMethod newton
            rng(1);
            K = objective(1e-3, 1e-3);
            [f0, g0] = K([0; 0]);
            opts.maxMethod = 'newton';
            [alpha, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, opts);
            testCase.verifyGreaterThan(fNew, f0);
            testCase.verifyLessThanOrEqual(abs(alpha - trueMaximum()), 0.05);
            testCase.verifyEqual(numel(unique(evals.A)), numel(evals.A));
        end

        function unknownMaxMethodGivesError(testCase)
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            opts.maxMethod = 'fzero';
            testCase.verifyError(@() line_search_main([0; 0], f0, g0, [2; 0.5], K, [], [], opts), ...
                ?MException);
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
            % förslaget efter de fyra förvalda punkterna blir 1 igen och
            % sökningen stoppar direkt.
            rng(1);
            K = @(x) deal(x(1) + 1e-3*randn, [1; 0] + 1e-3*randn(2, 1));
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K);
            testCase.verifyEqual(alpha, 1);
            testCase.verifyEqual(numel(evals.A), 5);
        end

        function initialPointsBeforeFirstFit(testCase)
            % alpha1, 3*alpha1/4, alpha1/2 och alpha1/4 mäts först, i den
            % ordningen, och inget polynom anpassas förrän alla fyra är
            % mätta (Q = NaN för de tre första). alpha1 = 100/S0 = 0.05 här.
            K = @(x) deal(-1e3*(x(1) - 1)^2, [-2e3*(x(1) - 1); 0]);
            [f0, g0] = K([0; 0]);              % S0 = 2000
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                struct('alphaMin', 1e-3));
            testCase.verifyEqual(evals.A(2:5), 0.05*[1 3/4 1/2 1/4], 'AbsTol', 1e-15);
            testCase.verifyTrue(all(isnan(evals.Q(1:4))));
            testCase.verifyFalse(isnan(evals.Q(5)));
            testCase.verifyEqual(evals.AHat(2:4), 0.05*[3/4 1/2 1/4], 'AbsTol', 1e-15);
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

        function maxLineSearchLimitsSimpleEvaluations(testCase)
            rng(1);
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            [~, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [2; 0.5], K, ...
                struct('maxLineSearch', 5), true);
            testCase.verifyEqual(numel(evals.A), 5);
            testCase.verifyTrue(contains(evals.reason, 'numel(alphas) = 5'), evals.reason);
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
            % Sökningen ska hoppa över NaN-punkterna, inte krascha, och returnera en ändlig
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

        function measuredWorseProposalStops(testCase)
            % F ökar ända fram till gränsen alpha = 1 (log(1 - alpha) = -Inf
            % där). Efter de förvalda punkterna föreslår polynomet gränsen
            % igen, en redan mätt ogiltig punkt. Sökningen ska stoppa och
            % returnera den bästa punkten, inte köra gränsen igen.
            Kb = @(x) deal(10*x(1) + 0.01*log(1 - x(1)), [10 - 0.01/(1 - x(1)); 0]);
            [f0, g0] = Kb([0; 0]);
            [alpha, ~, fNew, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], Kb, ...
                struct('alphaMin', 0.1));
            testCase.verifyTrue(contains(evals.reason, 'redan utvärderad'), evals.reason);
            testCase.verifyGreaterThan(fNew, f0);
            testCase.verifyEqual(fNew, max(evals.F));
            testCase.verifyLessThan(alpha, 1);
        end

        function neverStepsPastOne(testCase)
            % F ökar långt bortom alpha = 1, men steget får aldrig bli längre.
            K = @(x) deal(-0.01*(x(1) - 30)^2, [-0.02*(x(1) - 30); 0]);
            [f0, g0] = K([0; 0]);
            [alpha, ~, ~, ~, evals] = runQuiet([0; 0], f0, g0, [1; 0], K, ...
                struct('alphaMin', 0.1));
            testCase.verifyLessThanOrEqual(max(evals.A), 1);
            testCase.verifyEqual(alpha, 1);
        end

        function oldOptionsGiveError(testCase)
            % Parametrar från den gamla linjesökningen ska inte ignoreras tyst.
            K = objective(0.05, 0.1);
            [f0, g0] = K([0; 0]);
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, struct('minPoints', 5)), ?MException);
            % alphaMax och xLbd finns inte längre: taket är 1 och d är redan
            % skalad mot variansgränserna.
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, struct('alphaMax', 5)), ?MException);
            testCase.verifyError(@() runQuiet([0; 0], f0, g0, [2; 0.5], K, struct('xLbd', [0; 0])), ?MException);
        end


        function searchDirectionAndLineSearchRespectBounds(testCase)
            % Syntetisk, deterministisk målfunktion
            K = @(x) deal( ...
                -(x(1) + 1)^2 - x(2)^2, ...
                [-2*(x(1) + 1); -2*x(2)]);

            % Startpunkt nära undre gränsen
            x = [0.1; 0];
            xLbd = [0; -inf];

            [f, g] = K(x);

            % Steg 1: Beräkna sökriktning
            d = find_search_direction(x, g, xLbd, []);

            % Steg 2: Utför linjesökning
            [alpha, xNew, fNew] = runQuiet(x, f, g, d, K);

            % Den nya punkten ska respektera gränserna
            testCase.verifyGreaterThanOrEqual(xNew, xLbd);

            % Målfunktionen ska förbättras
            testCase.verifyGreaterThan(fNew, f);

            % Steglängden ska ligga inom tillåtet intervall
            testCase.verifyGreaterThan(alpha, 0);
            testCase.verifyLessThanOrEqual(alpha, 1);

            % De två funktionernas resultat ska vara konsistenta
            testCase.verifyEqual(xNew, x + alpha*d, ...
                'AbsTol', 1e-12);
        end

        function shrinkAlphaMaxCanBeDisabled(testCase)
            % Kvadratisk funktion med maximum vid alpha = 0.3.
            K = @(x) deal( ...
                -(x(1) - 0.3)^2, ...
                [-2*(x(1) - 0.3); 0]);

            x = [0; 0];
            d = [1; 0];
            [f, g] = K(x);

            % Kör med krympning aktiverad
            [~, ~, fShrink, ~, evalsShrink] = runQuiet( ...
                x, f, g, d, K, ...
                struct('shrinkAlphaMax', 1));

            % Kör med krympning avaktiverad
            [~, ~, fFixed, ~, evalsFixed] = runQuiet( ...
                x, f, g, d, K, ...
                struct('shrinkAlphaMax', 0));

            % Med krympning ska intervallet minska.
            testCase.verifyTrue( ...
                any([evalsShrink.fitAlphaMax{:}] < 1));

            % Utan krympning ska alphaMax alltid vara 1.
            testCase.verifyEqual( ...
                [evalsFixed.fitAlphaMax{:}], ...
                ones(size([evalsFixed.fitAlphaMax{:}])));

            % Båda metoderna ska hitta en förbättring.
            testCase.verifyGreaterThan(fShrink, f);
            testCase.verifyGreaterThan(fFixed, f);

            % Alla utvärderade steglängder ska ligga inom [0, 1].
            testCase.verifyTrue(all(evalsShrink.A >= 0 & evalsShrink.A <= 1));
            testCase.verifyTrue(all(evalsFixed.A >= 0 & evalsFixed.A <= 1));
        end

        function usesPointsOutsideShrunkAlphaMax(testCase)
            % Kvadratisk funktion med maximum vid alpha = 0.3.
            K = @(x) deal( ...
                -(x(1) - 0.3)^2, ...
                [-2*(x(1) - 0.3); 0]);

            x = [0; 0];
            d = [1; 0];
            [f, g] = K(x);

            % Aktivera krympning av alphaMax.
            [~, ~, ~, ~, evals] = runQuiet( ...
                x, f, g, d, K, ...
                struct('shrinkAlphaMax', 1, 'alphaMin', 1));

            % Kontrollera den första polynomanpassningen.
            alphaMax = evals.fitAlphaMax{1};
            t = evals.fitAlpha{1};

            % Sökintervallet ska ha krympt.
            testCase.verifyLessThan(alphaMax, 1);

            % Minst en mätpunkt ska ligga utanför
            % det nya sökintervallet.
            testCase.verifyTrue(any(t > 1));

            % Alla fem initiala mätpunkter ska användas
            % i den första polynomanpassningen.
            testCase.verifyEqual(numel(t), 5);

            % Kontrollera att även punkterna utanför
            % alphaMax finns med med rätt koordinater.
            testCase.verifyEqual( ...
                t * alphaMax, evals.A(1:5), ...
                'AbsTol', 1e-12);

            % Kontrollera att funktionsvärdena
            % från alla fem punkter används.
            testCase.verifyEqual( ...
                evals.fitF{1}, evals.F(1:5), ...
                'AbsTol', 1e-12);
        end

        function denoiseDirectionWorksWithLineSearch(testCase)
            % Konkav kvadratisk målfunktion med känt maximum.
            K = @(z) deal( ...
                -(z(1) - 0.8)^2 - (z(2) - 0.5)^2, ...
                [-2*(z(1) - 0.8); -2*(z(2) - 0.5)]);

            x = [0; 0];
            xLbd = -inf(2, 1);
            [f, g] = K(x);

            % 1. Beräkna vanlig sökriktning.
            [d, ~] = find_search_direction(x, g, xLbd, []);

            testCase.verifyGreaterThan(g' * d, 0);

            % 2. Skapa syntetisk mäthistorik med punkter
            % som spänner upp båda dimensionerna.
            X = [-0.5  0.5  0  0  0.4 -0.4  0.3  0;
                0    0   -0.5 0.5 -0.3  0.4  0.2  0];

            F = zeros(1, size(X, 2));
            G = zeros(2, size(X, 2));

            for i = 1:size(X, 2)
                [F(i), G(:, i)] = K(X(:, i));
            end

            evals.X = X;
            evals.F = F;
            evals.G = G;

            config.fixed_or_exact = 1;
            weights = struct('w_xi', 1, 'wbar_xi', 1);

            % 3. Beräkna denoisad sökriktning.
            dPrime = denoise(x, f, g, eye(2), ...
                {evals}, 2, [], [], [], xLbd, ...
                config, weights);

            testCase.verifyTrue(all(isfinite(dPrime)));
            testCase.verifyGreaterThan(norm(dPrime), 0);

            % Den nya riktningen ska vara en uppåtriktning.
            testCase.verifyGreaterThan(g' * dPrime, 0);

            % 4. Kör linjesökning med denoisad riktning.
            [alpha, xNew, fNew, gNew] = runQuiet( ...
                x, f, g, dPrime, K);

            % Linjesökningen ska hitta en förbättring.
            testCase.verifyGreaterThan(alpha, 0);
            testCase.verifyGreaterThan(fNew, f);
            testCase.verifyTrue(all(isfinite(gNew)));

            % Den nya punkten ska motsvara steget.
            testCase.verifyEqual(xNew, x + alpha*dPrime, ...
                'AbsTol', 1e-12);
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

function [alpha, xNew, fNew, gNew, evals] = runQuiet(x, f, g, d, K, opts, isSimple)
% line_search_main utan utskrifter.
% Om inget annat anges är alphaMin = 1, så första steget blir alpha = 1 och
% testerna ovan inte beror på initial_step_length.
if nargin < 6
    opts = struct();
end
if nargin < 7
    isSimple = 0;
end
if ~isfield(opts, 'alphaMin'), opts.alphaMin = 1; end
alpha = [];  xNew = [];  fNew = [];  gNew = [];  evals = [];
evalc('[alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, [], [], opts, isSimple);');
end
