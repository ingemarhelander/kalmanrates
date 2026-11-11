function tests = test_denoise_quadratic_approximation
% Tests the quadratic approximation used inside denoise using small examples.
% The approximation is a local function in denoise.m, so these tests call
% denoise with identity directions and check its estimated gradient.
    tests = functiontests(localfunctions);
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

function checkDenoising(testCase, f, g, seed)
% Sample a small cloud around one point and add modest measurement noise.
    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));
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
        eye(2), {evals}, 2, [], [], [], -inf(2, 1));

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
