classdef testFindSearchDirectionOrthogonal < matlab.unittest.TestCase

    methods (TestClassSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'AI_Implementations')));
        end
    end

    methods (Test)

        function initialDirectionIsNormalizedGradient(testCase)
            x = [0; 0];
            g = [3; 4];
            xLbd = [-inf; -inf];
            D = [];

            d = find_search_direction_orthogonal(x, g, xLbd, D);

            testCase.verifyEqual(d, [0.6; 0.8], 'AbsTol', 1e-12);
        end

        function directionIsOrthogonalToDependentNonorthogonalHistory(testCase)
            x = zeros(3, 1);
            g = [2; 3; 4];
            xLbd = -inf(3, 1);
            D = [1, 0, 1; 0, 1, 1; 1, 1, 2];

            d = find_search_direction_orthogonal(x, g, xLbd, D);

            testCase.verifyEqual(D' * d, zeros(size(D, 2), 1), ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(norm(d), 1, 'AbsTol', 1e-12);
            testCase.verifyGreaterThan(g' * d, 0);
        end

        function directionRespectsLowerBounds(testCase)
            x = [0.1; 2];
            g = [-2; 1];
            xLbd = [0; -inf];
            D = [];

            d = find_search_direction_orthogonal(x, g, xLbd, D);

            testCase.verifyGreaterThanOrEqual(x + d, xLbd);
            testCase.verifyGreaterThan(g' * d, 0);
        end

        function gradientInHistorySpanUsesGradientAndClearsHistory(testCase)
            x = zeros(3, 1);
            g = [3; 4; 0];
            xLbd = -inf(3, 1);
            D = [0.6; 0.8; 0];

            [d, Dnew] = find_search_direction_orthogonal(x, g, xLbd, D);

            testCase.verifyEqual(d, g / norm(g), 'AbsTol', 1e-12);
            testCase.verifyEmpty(Dnew);
        end
    end
end
