function tests = test_denoise_quadratic_approximation
% Tests the quadratic approximation used inside denoise using small examples.
% The approximation is a local function in denoise.m, so these tests call
% denoise with identity directions and check its estimated gradient.
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

function testQuadraticFunction(testCase)
    [f, g] = makeFunction(1);
    checkDenoising(testCase, f, g, 11);
end

function testSmoothCubicFunction(testCase)
    [f, g] = makeFunction(2);
    checkDenoising(testCase, f, g, 22);
end

function testSmoothSineFunction(testCase)
    [f, g] = makeFunction(3);
    checkDenoising(testCase, f, g, 33);
end

function testFixedPointMatchesDirect(testCase)
    % Exakt kvadratisk funktion utan brus
    [f, g] = makeFunction(1);

    center = [0.35; -0.2];

    X = [0 1 -1 0 0.5 -0.5 0.3;
         0 0  0 1 -1   0.5 -0.2];

    F = zeros(1, size(X, 2));
    G = zeros(2, size(X, 2));

    for i = 1:size(X, 2)
        F(i) = f(X(:, i));
        G(:, i) = g(X(:, i));
    end

    % Deterministiskt brus för reproducerbart test
    F = F + 0.001 * sin(1:numel(F));
    G = G + 0.001 * [sin(1:size(G,2)); cos(1:size(G,2))];

    evals.X = X;
    evals.F = F;
    evals.G = G;

    weights = struct('w_xi', 1, 'wbar_xi', 1);
    xLbd = -inf(2, 1);

    % Fixed-point iteration
    config.fixed_or_exact = 0;
    dFixed = denoise(center, f(center), g(center), ...
        eye(2), {evals}, 2, [], [], [], ...
        xLbd, config, weights);

    % Direkt lösning
    config.fixed_or_exact = 1;
    dDirect = denoise(center, f(center), g(center), ...
        eye(2), {evals}, 2, [], [], [], ...
        xLbd, config, weights);

    % Metoderna ska ge samma sökriktning
    verifyEqual(testCase, dFixed, dDirect, ...
        'AbsTol', 1e-6);
end

function testInsufficientSamples(testCase)
    % Endast en mätpunkt i två dimensioner
    x = [0.5; 0.3];
    f = 2;
    g = [1; -1];

    evals.X = x;
    evals.F = f;
    evals.G = g;

    B = eye(2);
    xLbd = -inf(2, 1);

    config.fixed_or_exact = 1;
    weights = struct('w_xi', 1, 'wbar_xi', 1);

    d = denoise(x, f, g, B, {evals}, 2, ...
        [], [], [], xLbd, config, weights);

    % Ingen unik kvadratisk approximation är möjlig
    verifyEqual(testCase, d, zeros(2, 1));
end

function testDenoiseRespectsLowerBounds(testCase)
    % Startpunkt nära undre gränsen
    x = [0.1; 0];
    xLbd = [0; -inf];

    % Syntetiska mätpunkter, med startpunkten sist
    X = [0 1 -1 0 0.5 -0.5 0.3 0.1;
         0 0  0 1 -1   0.5 -0.2 0];

    % Kvadratisk målfunktion
    f = @(z) -(z(1) + 1)^2 - z(2)^2;
    g = @(z) [-2*(z(1) + 1); -2*z(2)];

    F = zeros(1, size(X, 2));
    G = zeros(2, size(X, 2));

    for i = 1:size(X, 2)
        F(i) = f(X(:, i));
        G(:, i) = g(X(:, i));
    end

    evals.X = X;
    evals.F = F;
    evals.G = G;

    config.fixed_or_exact = 1;
    weights = struct('w_xi', 1, 'wbar_xi', 1);

    d = denoise(x, f(x), g(x), eye(2), {evals}, 2, ...
        [], [], [], xLbd, config, weights);

    % Riktningen ska vara ändlig och inte noll
    verifyTrue(testCase, all(isfinite(d)));
    verifyGreaterThan(testCase, norm(d), 0);

    % Riktningen ska peka mot den undre gränsen
    verifyLessThan(testCase, d(1), 0);

    % Ett helt steg får inte överskrida gränsen
    verifyGreaterThanOrEqual(testCase, x + d, xLbd);
end

function testCollinearSamples(testCase)
    % Flera mätpunkter, men alla ligger på samma linje
    X = [-2 -1 0 1 2;
          0  0 0 0 0];

    f = @(x) -(x(1)-1)^2 - x(2)^2;
    g = @(x) [-2*(x(1)-1); -2*x(2)];

    F = zeros(1, size(X, 2));
    G = zeros(2, size(X, 2));

    for i = 1:size(X, 2)
        F(i) = f(X(:, i));
        G(:, i) = g(X(:, i));
    end

    evals.X = X;
    evals.F = F;
    evals.G = G;

    x = X(:, end);
    config.fixed_or_exact = 1;
    weights = struct('w_xi', 1, 'wbar_xi', 1);

    d = denoise(x, F(end), G(:, end), ...
        eye(2), {evals}, 2, [], [], [], ...
        -inf(2,1), config, weights);

    % Punkterna täcker inte båda dimensionerna
    verifyEqual(testCase, d, zeros(2,1));
end

function testMultipleEvaluationHistories(testCase)
    % Kvadratisk målfunktion
    [f, g] = makeFunction(1);

    X = [0 1 -1 0 0.5 -0.5 0.3 0.1;
         0 0  0 1 -1  0.5 -0.2 0.2];

    F = zeros(1, size(X, 2));
    G = zeros(2, size(X, 2));

    for i = 1:size(X, 2)
        F(i) = f(X(:, i));
        G(:, i) = g(X(:, i));
    end

    % Alla observationer i en historikpost
    evalsAll.X = X;
    evalsAll.F = F;
    evalsAll.G = G;

    % Samma observationer fördelade på två historikposter
    evals1.X = X(:, 1:4);
    evals1.F = F(1:4);
    evals1.G = G(:, 1:4);

    evals2.X = X(:, 5:8);
    evals2.F = F(5:8);
    evals2.G = G(:, 5:8);

    x = X(:, end);
    config.fixed_or_exact = 1;
    weights = struct('w_xi', 1, 'wbar_xi', 1);

    % Kör med en respektive två historikposter
    dAll = denoise(x, F(end), G(:,end), ...
        eye(2), {evalsAll}, 2, [], [], [], ...
        -inf(2,1), config, weights);

    dSplit = denoise(x, F(end), G(:,end), ...
        eye(2), {evals1, evals2}, 2, [], [], [], ...
        -inf(2,1), config, weights);

    % Samma mätpunkter ska ge samma resultat
    verifyEqual(testCase, dSplit, dAll, ...
        'AbsTol', 1e-10);

    % Resultatet ska vara en användbar sökriktning
    verifyTrue(testCase, all(isfinite(dSplit)));
    verifyEqual(testCase, norm(dSplit), 1, ...
        'AbsTol', 1e-10);
end

function checkDenoising(testCase, f, g, seed)
% Sample a small cloud around one point and add modest measurement noise.
    previousRng = rng;
    rngCleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(seed, 'twister');

    center = [0.35; -0.2];
    sampleCount = 31;
    cloudRadius = 0.12;

    % Keep the centre as the final sample because the model is centred there.
    X = [center + cloudRadius * randn(2, sampleCount - 1), center];
    F = zeros(1, sampleCount);
    G = zeros(2, sampleCount);

    for sample = 1:sampleCount
        F(sample) = f(X(:, sample));
        G(:, sample) = g(X(:, sample));
    end

    % Add independent noise to the function and gradient measurements.
    trueGradient = g(center);
    functionNoise = 0.01 * randn(size(F));
    gradientNoise = 0.08 * randn(size(G));

    % Make the centre gradient noticeably noisy so the test can check that
    % information from the surrounding cloud improves its direction.
    perpendicular = [-trueGradient(2); trueGradient(1)];
    perpendicular = perpendicular / norm(perpendicular);
    gradientNoise(:, end) = 0.6 * norm(trueGradient) * perpendicular;

    noisyF = F + functionNoise;
    noisyG = G + gradientNoise;

    % denoise's other inputs are not needed for this isolated test. Identity
    % directions make it use the full two-dimensional sample cloud.
    evals.X = X;
    evals.F = noisyF;
    evals.G = noisyG;
    estimatedDirection = denoise(center, noisyF(end), noisyG(:, end), ...
        eye(2), {evals}, 2, [], [], [], -inf(2, 1), ...
        struct('fixed_or_exact', 1), struct('w_xi', 1, 'wbar_xi', 1));

    trueDirection = trueGradient / norm(trueGradient);
    noisyCenterDirection = noisyG(:, end) / norm(noisyG(:, end));

    % Compare directions, not gradient magnitudes: denoise returns a unit
    % search direction rather than the fitted gradient itself.
    rawError = norm(noisyCenterDirection - trueDirection);
    denoisedError = norm(estimatedDirection - trueDirection);

    verifyLessThan(testCase, denoisedError, rawError);
    verifyLessThan(testCase, denoisedError, 0.25);
end

function [f, g] = makeFunction(whichFunction)
% Three inexpensive objectives: quadratic and two smooth nonquadratics.
    H = [2.0, 0.3; 0.3, 1.4];
    b = [0.8; 0.6];

    switch whichFunction
        case 1
            f = @(x) 0.5 * x' * H * x + b' * x;
            g = @(x) H * x + b;
        case 2
            f = @(x) 0.5 * x' * H * x + b' * x + ...
                0.1 * (x(1)^3 - 0.5 * x(2)^3);
            g = @(x) H * x + b + [0.3 * x(1)^2; -0.15 * x(2)^2];
        case 3
            f = @(x) 0.5 * x' * H * x + b' * x + ...
                0.15 * sin(x(1)) + 0.1 * cos(x(2));
            g = @(x) H * x + b + [0.15 * cos(x(1)); -0.1 * sin(x(2))];
        otherwise
            error('Unknown test function number.');
    end
end
