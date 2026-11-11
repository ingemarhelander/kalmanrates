classdef testSelectPolyDegreeLr < matlab.unittest.TestCase
% testSelectPolyDegreeLr - Tester för gradval med LR-test
% (select_poly_degree_lr, lr_statistic och ml_poly_main med 'lr').
%
% Kör från projektroten med:  setupPaths; runtests('tests/testSelectPolyDegreeLr.m')

    methods (TestClassSetup)
        function addPaths(testCase)
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
        function quadraticStaysTwo(testCase)
            % Rätt grad och bara brus: graden ska sällan höjas.
            q = selectMany(@(a) 1 + 3*a - 5*a.^2, @(a) 3 - 10*a, 8);
            testCase.verifyGreaterThanOrEqual(mean(q == 2), 0.9);
        end

        function cubicGivesThree(testCase)
            q = selectMany(@(a) 1 + 2*a - 6*a.^3, @(a) 2 - 18*a.^2, 6);
            testCase.verifyGreaterThanOrEqual(mean(q == 3), 0.9);
        end

        function quarticGivesFour(testCase)
            q = selectMany(@(a) 1 - 4*a.^3 + 5*a.^4, @(a) -12*a.^2 + 20*a.^3, 6);
            testCase.verifyGreaterThanOrEqual(mean(q == 4), 0.9);
        end

        function maxDegreeIsMMinusTwo(testCase)
            % Fjärdegrads med för få punkter: graden stannar vid m - 2 (minst 2).
            for m = 2:5
                q = selectMany(@(a) 1 - 4*a.^3 + 5*a.^4, @(a) -12*a.^2 + 20*a.^3, m);
                testCase.verifyTrue(all(q <= max(2, m - 2)), sprintf('m = %d', m));
            end
        end

        function exactPolynomialGivesDegreeTwo(testCase)
            alpha = linspace(0, 1, 8)';
            q = [];
            evalc('q = select_poly_degree_lr(alpha, 1 + 3*alpha - 5*alpha.^2, 3 - 10*alpha);');
            testCase.verifyEqual(q, 2);
        end

        function mlPolyMainChoosesMethod(testCase)
            % ml_poly_main med 'lr' ger samma som select_poly_degree_lr,
            % och en okänd metod ger fel.
            rng(1);
            alpha = linspace(0, 1, 8)';
            F = 1 - 4*alpha.^3 + 5*alpha.^4 + 0.003*randn(8, 1);
            S = -12*alpha.^2 + 20*alpha.^3 + 0.015*randn(8, 1);
            q1 = [];  q2 = [];
            evalc('q1 = ml_poly_main(alpha, F, S, [], ''lr'');');
            evalc('q2 = select_poly_degree_lr(alpha, F, S);');
            testCase.verifyEqual(q1, q2);
            testCase.verifyError(@() ml_poly_main(alpha, F, S, [], 'xyz'), ?MException);
        end
    end
end

function q = selectMany(F_true, S_true, m)
% Vald grad för 50 slumpade dataset med m jämnt fördelade punkter.
    alpha = linspace(0, 1, m)';
    q = zeros(50, 1);
    for s = 1:50
        rng(s);
        F = F_true(alpha) + 0.003*randn(m, 1);
        S = S_true(alpha) + 0.015*randn(m, 1);
        qs = [];
        evalc('qs = select_poly_degree_lr(alpha, F, S);');
        q(s) = qs;
    end
end
