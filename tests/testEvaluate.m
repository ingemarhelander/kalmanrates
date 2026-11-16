classdef testEvaluate < matlab.unittest.TestCase
% testEvaluate - Tester för avbrottskriterierna i linjesökningen (evaluate,
% komponent 10).
%
% evaluate anropas direkt med ett handskrivet polynom, så varken Kalman eller
% ML-anpassningen behövs. Varje test utgår från baseline(), där inget
% kriterium är uppfyllt, och ändrar en sak i taget.
%
% Kör från projektroten med:  runtests('tests/testEvaluate.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(testDir), 'scripts', 'line_search')));
        end
    end

    methods (Test)
        % ------------------------------ Grundfall ------------------------------
        function baselineDoesNotStop(testCase)
            [stop, reason] = run(baseline());
            testCase.verifyFalse(stop);
            testCase.verifyEmpty(reason);
        end

        % ------------------------- Kriterium 1, lutning ------------------------
        function flatSlopeStops(testCase)
            a = baseline();
            a.fj = 0.1;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'lutningen');
        end

        function flatSlopeWithoutImprovementDoesNotStop(testCase)
            a = baseline();
            a.fj = 0.1;
            a.F = [0, -0.1];
            testCase.verifyFalse(run(a));
        end

        function slopeEqualToLimitStops(testCase)
            % Villkoret är <=, så likhet stoppar.
            a = baseline();
            a.fj = a.opts.c*a.f0;
            testCase.verifyTrue(run(a));
        end

        function slopeSignDoesNotMatter(testCase)
            a = baseline();
            a.fj = -0.1;
            testCase.verifyTrue(run(a), 'negativ lutning i alpha_j');
            a.fj = -0.8;
            testCase.verifyFalse(run(a), 'brant negativ lutning i alpha_j');
            a.fj = 0.1;
            a.f0 = -1;
            testCase.verifyTrue(run(a), 'negativ lutning i startpunkten');
        end

        % ------------------- Kriterium 2, för nära (Req. 31) -------------------
        function tooCloseStops(testCase)
            a = baseline();
            a.alphaHat = 0.2 + 1e-5;       % gränsen är 1e-3*alphaHat = 2e-4
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'Req. 31');
        end

        function tooCloseWithoutImprovementDoesNotStop(testCase)
            % Det fallet tar line_search_main hand om (kortare steg).
            a = baseline();
            a.alphaHat = 0.2 + 1e-5;
            a.F = [0, -0.1];
            testCase.verifyFalse(run(a));
        end

        function tooCloseToInvalidPointStops(testCase)
            % En punkt som gett NaN räknas också som utvärderad.
            a = baseline();
            a.alphas = [0, 0.4, 0.2];
            a.F = [0, NaN, 0.18];
            a.ok = [true, false, true];
            a.alphaHat = 0.4 + 1e-5;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'Req. 31');
        end

        function distanceIsRelativeToProposedStep(testCase)
            % Största steget 1 var ett misslyckat långt steg, förbättringen
            % finns vid 0.005. Avståndet 5e-4 dit är inte för nära (gränsen
            % är 1e-3*0.0055 = 5.5e-6, inte 1e-3*1), men avståndet 1e-6 är det.
            a = baseline();
            a.opts.noiseFactor = 0;        % stäng av kriterium 3
            a.alphas = [0, 1, 0.005];
            a.F = [0, -13, 0.05];
            a.ok = [true, true, true];
            a.alphaHat = 0.0055;
            testCase.verifyFalse(run(a));
            a.alphaHat = 0.005 + 1e-6;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'för nära');
        end

        function tooCloseToStartPointStops(testCase)
            % alphaHat under dMinRel gånger det minsta utvärderade steget
            % räknas som startpunkten (gränsen är 1e-3*1e-3 = 1e-6), så att
            % sökningen inte kryper mot alpha = 0.
            a = baseline();
            a.opts.noiseFactor = 0;        % stäng av kriterium 3
            a.alphas = [0, 1e-3];
            a.F = [0, 1e-3];
            a.alphaHat = 1e-5;
            testCase.verifyFalse(run(a));
            a.alphaHat = 1e-7;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'Req. 31');
        end

        function tooCloseToStartPointWithoutImprovementDoesNotStop(testCase)
            % Det fallet tar line_search_main hand om.
            a = baseline();
            a.opts.noiseFactor = 0;
            a.alphas = [0, 1e-3];
            a.F = [0, -1e-3];
            a.alphaHat = 1e-12;
            testCase.verifyFalse(run(a));
        end

        % -------------------------- Kriterium 3, brus --------------------------
        function smallPredictedIncreaseStops(testCase)
            % Förutsagd ökning p(1) - p(0.2) = 0.32.
            a = baseline();
            a.fit.sigmaF = 0.5;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'bruset');
        end

        function smallPredictedIncreaseWithoutImprovementDoesNotStop(testCase)
            a = baseline();
            a.fit.sigmaF = 0.5;
            a.F = [0, -0.1];
            testCase.verifyFalse(run(a));
        end

        function predictedDecreaseStops(testCase)
            % p(0.1) < p(0.2): polynomet förutsäger ingen ökning alls.
            a = baseline();
            a.alphaHat = 0.1;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'bruset');
        end

        function largerNoiseFactorStopsEarlier(testCase)
            % Förutsagd ökning 0.32, sigmaF = 0.2.
            a = baseline();
            a.fit.sigmaF = 0.2;
            a.opts.noiseFactor = 1;
            testCase.verifyFalse(run(a));
            a.opts.noiseFactor = 2;
            testCase.verifyTrue(run(a));
        end

        function exactFitDoesNotStopOnNoise(testCase)
            a = baseline();
            a.fit.sigmaF = 0;
            testCase.verifyFalse(run(a));
        end

        function noiseComparesWithValidPointsOnly(testCase)
            % Den ogiltiga punkten alpha = 1 är polynomets maximum. Räknades
            % den med skulle den förutsagda ökningen bli negativ och ge stopp.
            a = baseline();
            a.alphas = [0, 1, 0.2];
            a.F = [0, NaN, 0.18];
            a.ok = [true, false, true];
            a.alphaHat = 0.9;
            testCase.verifyFalse(run(a));
        end

        % ------------------------------ Gemensamt ------------------------------
        function minPointsBlocksSlopeAndNoise(testCase)
            a = baseline();
            a.opts.minPoints = 3;          % bara 2 punkter utvärderade
            b = a;  b.fj = 0.1;
            testCase.verifyFalse(run(b), 'lutning');
            b = a;  b.fit.sigmaF = 0.5;
            testCase.verifyFalse(run(b), 'brus');
        end

        function minPointsDoesNotBlockTooClose(testCase)
            a = baseline();
            a.opts.minPoints = 3;
            a.alphaHat = 0.2 + 1e-5;
            [stop, reason] = run(a);
            testCase.verifyTrue(stop);
            testCase.verifySubstring(reason, 'Req. 31');
        end

        function minPointsReachedAllowsStop(testCase)
            % Startpunkten räknas med: 3 punkter räcker för minPoints = 3.
            a = baseline();
            a.opts.minPoints = 3;
            a.alphas = [0, 0.1, 0.2];
            a.F = [0, 0.095, 0.18];
            a.ok = [true, true, true];
            a.fj = 0.1;
            testCase.verifyTrue(run(a));
        end

        function firstFulfilledCriterionGivesReason(testCase)
            % Ordningen är lutning, för nära, brus.
            a = baseline();
            a.fj = 0.1;
            a.alphaHat = 0.2 + 1e-5;
            a.fit.sigmaF = 0.5;
            [~, reason] = run(a);
            testCase.verifySubstring(reason, 'lutningen');
            a.fj = 0.8;
            [~, reason] = run(a);
            testCase.verifySubstring(reason, 'Req. 31');
        end

        function invalidPointIsNotAnImprovement(testCase)
            % F = Inf i en ogiltig punkt räknas inte som att F har ökat.
            a = baseline();
            a.alphas = [0, 0.4, 0.2];
            a.F = [0, Inf, -0.1];
            a.ok = [true, false, true];
            a.fj = 0.1;
            testCase.verifyFalse(run(a));
        end
    end
end

function a = baseline()
% Indata till evaluate där inget kriterium är uppfyllt.
% p(alpha) = alpha - 0.5*alpha^2, med maximum p(1) = 0.5. Två punkter är
% utvärderade, alpha = 0 och 0.2, och F har ökat. Nästa steg är alphaHat = 1.
%   1. lutning:  |f'_j| = 0.8 > c*|f'_0| = 0.2
%   2. för nära: avståndet 0.8 > dMinRel*alphaHat = 1e-3
%   3. brus:     förutsagd ökning p(1) - p(0.2) = 0.32 > sigmaF = 0.01
    a.fit.y = [0; 1; -0.5];
    a.fit.sigmaF = 0.01;
    a.alphaHat = 1;
    a.alphas = [0, 0.2];
    a.ok = [true, true];
    a.F = [0, 0.18];
    a.fj = 0.8;
    a.f0 = 1;
    a.opts = struct('c', 0.2, 'dMinRel', 1e-3, 'noiseFactor', 1, 'minPoints', 0);
end

function [stop, reason] = run(a)
    [stop, reason] = evaluate(a.fit, a.alphaHat, a.alphas, a.ok, a.F, a.fj, a.f0, a.opts);
end
