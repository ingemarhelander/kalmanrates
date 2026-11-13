function tests = testFindSearchDirection
% Tester för find_search_direction: ortogonalisering mot tidigare riktningar.
%
% Kör från projektroten med:  setupPaths; runtests('tests/testFindSearchDirection.m')
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.TestData.originalPath = path;
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(projectRoot, 'scripts', 'functions'));
end

function teardownOnce(testCase)
    path(testCase.TestData.originalPath);
end

function testOrthogonalToScaledDirections(testCase)
% Kolumnerna i D är skalade (som efter gränskontrollen) och inte inbördes
% ortogonala. Riktningen ska ändå bli ortogonal mot hela span(D).
    rng(4, 'twister');
    n = 6;
    x = ones(n, 1);
    xLbd = -inf(n, 1);              % inga gränser, så d skalas inte om
    g = 5e8 * randn(n, 1);          % stor gradient, som i Kalmanproblemet
    D = [1e-3 * randn(n, 1), 2e-4 * randn(n, 1), 7 * randn(n, 1)];
    D(:, 2) = D(:, 2) + 0.5 * D(:, 1);   % inte ortogonala

    d = runQuiet(x, g, xLbd, D);

    testCase.verifyEqual(norm(d), 1, 'AbsTol', 1e-12);
    testCase.verifyLessThan(max(abs(orth(D)' * d)), 1e-12);
    testCase.verifyGreaterThan(g' * d, 0);   % fortfarande en ascent-riktning
end

function testGradientInSpanFallsBackToUnitGradient(testCase)
% Händer efter ett misslyckat steg: x och g är oförändrade och g ligger i
% span(D), eftersom senaste riktningen är g minus dess projektion på de
% äldre riktningarna.
    n = 4;
    x = ones(n, 1);
    xLbd = -inf(n, 1);
    g = [3; -1; 0; 0];
    D = [1e-3 * [1; 0; 0; 0], 1e-3 * [0; 1; 0; 0]];   % g ligger i span(D)

    d = runQuiet(x, g, xLbd, D);

    testCase.verifyEqual(d, g / norm(g), 'AbsTol', 1e-12);
end

function testBoundScalingKeepsDirection(testCase)
% Med en nedre gräns skalas d så att alpha = 1 är tillåtet, men riktningen
% är fortfarande den ortogonaliserade gradienten.
    x = [1e-3; 1; 1];
    xLbd = [0; -inf; -inf];
    g = [-2; 1; 0.5];
    D = 1e-2 * [0; 0; 1];

    d = runQuiet(x, g, xLbd, D);

    dUnit = [-2; 1; 0] / norm([-2; 1; 0]);
    testCase.verifyEqual(d / norm(d), dUnit, 'AbsTol', 1e-12);
    testCase.verifyGreaterThanOrEqual(x(1) + d(1), 0);   % alpha = 1 är tillåtet
end

function d = runQuiet(x, g, xLbd, D)
    d = [];
    evalc('d = find_search_direction(x, g, xLbd, D);');
end
