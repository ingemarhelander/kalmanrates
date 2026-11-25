classdef testFindStepLengthNewton < matlab.unittest.TestCase
% testFindStepLengthNewton - Tester för find_step_length_newton och
% poly_curvature på kända polynom.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testFindStepLengthNewton.m')

    methods (TestClassSetup)
        function addPaths(testCase)
            testDir = fileparts(mfilename('fullpath'));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(fileparts(testDir), 'scripts', 'line_search')));
        end
    end

    methods (Test)
        function curvatureOfKnownPolynomial(testCase)
            % p = 1 + 2a + 3a^2 + 4a^3: p'' = 6 + 24a
            fit.y = [1; 2; 3; 4];
            a = [0, 0.5, 2];
            testCase.verifyEqual(poly_curvature(fit, a), 6 + 24*a, 'AbsTol', 1e-12);
            fit.y = [1; 2];                        % linjärt: p'' = 0
            testCase.verifyEqual(poly_curvature(fit, a), zeros(size(a)));
        end

        function findsQuadraticMaximumExactly(testCase)
            % Grid ger bara h/2 här, Newton ska ge maskinprecision
            fit.y = coeffs(0.3047);
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 0.3047, 'AbsTol', 1e-12);
        end

        function findsSmallMaximumExactly(testCase)
            % Maximum långt under första punkten 1/50 = 0.02. Grid med
            % h = 1e-3 skulle ge 5e-4.
            fit.y = coeffs(1e-4);
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 1e-4, 'AbsTol', 1e-12);
        end

        function takesFirstMaximumOfCubicSlope(testCase)
            % p' = -(a - 0.2)(a - 0.5)(a - 0.8): max vid 0.2 och 0.8
            p = polyint(-poly([0.2 0.5 0.8]));
            fit.y = flipud(p(:));
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 0.2, 'AbsTol', 1e-12);
        end

        function findsMaximumOfNonQuadraticPolynomial(testCase)
            % p' = (0.3 - a)(1 + 5a^2) > 0 före 0.3, < 0 efter
            dp = conv([-1 0.3], [5 0 1]);
            p = polyint(dp);
            fit.y = flipud(p(:));
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 0.3, 'AbsTol', 1e-12);
        end

        function noMaximumGivesAlphaMax(testCase)
            fit.y = [0; 1; 0.5];                   % p' > 0 för alpha > 0
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 1);
            testCase.verifyEqual(find_step_length_newton(fit, 50, 0.37), 0.37);
        end

        function fallingFromStartGivesHalfFirstStep(testCase)
            fit.y = [0; -1; -1];                   % p'(0) = -1
            testCase.verifyEqual(find_step_length_newton(fit, 50, 1), 0.01, 'AbsTol', 1e-15);
        end

        function agreesWithGridWithinHalfStep(testCase)
            % Samma maximum som find_step_length, upp till dess noggrannhet h/2
            h = 1e-3;
            for a = [0.013, 0.3047, 0.777]
                fit.y = coeffs(a);
                aGrid = find_step_length(fit, h, 1);
                aNewton = find_step_length_newton(fit, 50, 1);
                testCase.verifyLessThanOrEqual(abs(aGrid - aNewton), h/2 + 1e-12);
            end
        end
    end
end

function y = coeffs(a)
% Koefficienter (stigande potenser) för p = -(alpha - a)^2.
    y = [-a^2; 2*a; -1];
end
