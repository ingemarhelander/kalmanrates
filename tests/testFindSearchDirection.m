
classdef testFindSearchDirection < matlab.unittest.TestCase

    methods (TestClassSetup)
        function addPaths(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture( ...
                matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'scripts', 'functions')));
        end
    end

    methods (Test)

        function initialDirectionIsNormalizedGradient(testCase)
            x = [0; 0];
            g = [3; 4];
            xLbd = [-inf; -inf];
            D = [];

            d = find_search_direction(x, g, xLbd, D);

            testCase.verifyEqual(d, [0.6; 0.8], ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(norm(d), 1, ...
                'AbsTol', 1e-12);
        end

        function newDirectionIsOrthogonal(testCase)
            x = zeros(3, 1);
            g = [2; 3; 4];
            xLbd = -inf(3, 1);
            D = [1; 0; 0];

            d = find_search_direction(x, g, xLbd, D);

            testCase.verifyEqual(D' * d, 0, ...
                'AbsTol', 1e-12);
            testCase.verifyEqual(norm(d), 1, ...
                'AbsTol', 1e-12);
        end

        function directionRespectsLowerBounds(testCase)
            x = [0.1; 2];
            g = [-2; 1];
            xLbd = [0; -inf];
            D = [];

            d = find_search_direction(x, g, xLbd, D);

            % alpha = 1 ska vara tillåtet
            testCase.verifyGreaterThanOrEqual( ...
                x + d, xLbd);

            % Sökriktningen ska fortfarande vara en uppåtriktning
            testCase.verifyGreaterThan(g' * d, 0);
        end

        function gradientInPreviousSpanRequestsDenoising(testCase)
            x = zeros(3, 1);
            g = [3; 4; 0];
            xLbd = -inf(3, 1);

            % Gradientens riktning finns redan i D
            D = [0.6; 0.8; 0];

            [d, Dnew, needsDenoising] = find_search_direction(x, g, xLbd, D);

            % Algoritmen ska återgå till normaliserad gradient
            testCase.verifyEqual(d, g / norm(g), ...
                'AbsTol', 1e-12);

            testCase.verifyEqual(Dnew, D);
            testCase.verifyTrue(needsDenoising);

            % Riktningen ska fortfarande vara en uppåtriktning
            testCase.verifyGreaterThan(g' * d, 0);
        end
    end
end