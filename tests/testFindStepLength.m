classdef testFindStepLength < matlab.unittest.TestCase
% testFindStepLength - Tester för find_step_length på kända polynom.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testFindStepLength.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(testDir), 'scripts', 'line_search')));
        end
    end

    methods (Test)
        function findsMaximumWithinHalfStep(testCase)
            % p = -(alpha - 0.3047)^2: maximum vid 0.3047
            fit.y = coeffs(0.3047);
            h = 0.01;
            alpha = find_step_length(fit, h, 1);
            testCase.verifyLessThanOrEqual(abs(alpha - 0.3047), h/2 + 1e-12);
        end

        function returnsMidpointOfTheTwoLastSteps(testCase)
            % Maximum vid 0.305: p' blir negativ i 0.31, föregående steg 0.30
            fit.y = coeffs(0.305);
            alpha = find_step_length(fit, 0.01, 1);
            testCase.verifyEqual(alpha, 0.305, 'AbsTol', 1e-12);
        end

        function takesFirstMaximum(testCase)
            % p'(alpha) = (alpha - 0.2)(alpha - 0.5)(alpha - 0.8)*(-1):
            % maximum vid 0.2 och 0.8, minimum vid 0.5. Det första ska väljas.
            r = [0.2 0.5 0.8];
            dp = -poly(r);                         % p' med fallande potenser
            p = polyint(dp);                       % p med fallande potenser
            fit.y = flipud(p(:));                  % stigande potenser som i poly_fit_ml
            alpha = find_step_length(fit, 1e-3, 1);
            testCase.verifyEqual(alpha, 0.2, 'AbsTol', 1e-3);
        end

        function noMaximumGivesAlphaMax(testCase)
            fit.y = [0; 1; 0.5];                   % p' > 0 för alpha > 0
            testCase.verifyEqual(find_step_length(fit, 0.01, 1), 1);
            testCase.verifyEqual(find_step_length(fit, 0.01, 0.37), 0.37);
        end

        function maximumBeforeFirstStepGivesHalfStep(testCase)
            fit.y = coeffs(0.001);
            testCase.verifyEqual(find_step_length(fit, 0.01, 1), 0.005, 'AbsTol', 1e-12);
        end
    end
end

function y = coeffs(a)
% Koefficienter (stigande potenser) för p = -(alpha - a)^2.
    y = [-a^2; 2*a; -1];
end
