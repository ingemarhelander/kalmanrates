function results = demo_denoise
% DEMO_DENOISE Visualize denoising on a simple noisy quadratic.
% Run from the project root with: demo_denoise
% The figure compares one sample cloud with angular errors from repeated
% noise realizations. A smaller angle means a direction closer to truth.

    originalPath = path;
    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(projectRoot, 'scripts', 'functions'));
    pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
    previousRng = rng;
    rngCleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
    rng(42, 'twister');

    H = [2.0, 0.3; 0.3, 1.4];
    b = [0.8; 0.6];
    f = @(x) 0.5 * x' * H * x + b' * x;
    g = @(x) H * x + b;
    center = [0.35; -0.2];
    trueDirection = g(center) / norm(g(center));

    sampleCount = 31;
    trialCount = 100;
    cloudRadius = 0.12;
    functionNoiseStd = 0.01;
    gradientNoiseStd = 0.08;
    X = [center + cloudRadius * randn(2, sampleCount - 1), center];

    F = zeros(1, sampleCount);
    G = zeros(2, sampleCount);
    for sample = 1:sampleCount
        F(sample) = f(X(:, sample));
        G(:, sample) = g(X(:, sample));
    end

    perpendicular = [-trueDirection(2); trueDirection(1)];
    rawErrors = zeros(trialCount, 1);
    denoisedErrors = zeros(trialCount, 1);
    representative = struct();

    for trial = 1:trialCount
        noisyF = F + functionNoiseStd * randn(size(F));
        noisyG = G + gradientNoiseStd * randn(size(G));
        noisyG(:, end) = G(:, end) + ...
            0.6 * norm(G(:, end)) * perpendicular;

        evals = struct('X', X, 'F', noisyF, 'G', noisyG);
        config = struct('fixed_or_exact', 1);
        weights = struct('w_xi', 1, 'wbar_xi', 1);

        % denoise reports model conditioning; capture that text for this demo.
        evalc(['denoised = denoise(center, noisyF(end), noisyG(:, end), ' ...
            'eye(2), {evals}, 2, [], [], [], -inf(2, 1), config, weights);']);

        rawDirection = noisyG(:, end) / norm(noisyG(:, end));
        rawErrors(trial) = directionErrorDegrees(rawDirection, trueDirection);
        denoisedErrors(trial) = directionErrorDegrees(denoised, trueDirection);

        if trial == 1
            representative.noisyG = noisyG(:, end);
            representative.denoised = denoised;
        end
    end

    improvedPercent = 100 * mean(denoisedErrors < rawErrors);
    fprintf('Noisy quadratic demo (%d trials)\n', trialCount);
    fprintf('Mean angular error:   raw %5.1f deg, denoised %5.1f deg\n', ...
        mean(rawErrors), mean(denoisedErrors));
    fprintf('Median angular error: raw %5.1f deg, denoised %5.1f deg\n', ...
        median(rawErrors), median(denoisedErrors));
    fprintf('Denoising improved the direction in %.0f%% of trials.\n', ...
        improvedPercent);

    results = struct('rawErrors', rawErrors, ...
        'denoisedErrors', denoisedErrors, 'improvedPercent', improvedPercent);
    plotDenoiseResults(H, b, center, X, trueDirection, representative, ...
        rawErrors, denoisedErrors);
end

function errorDegrees = directionErrorDegrees(direction, trueDirection)
    cosine = max(-1, min(1, direction' * trueDirection));
    errorDegrees = acosd(cosine);
end

function plotDenoiseResults(H, b, center, X, trueDirection, representative, ...
        rawErrors, denoisedErrors)
    figure('Name', 'Denoising on a noisy quadratic', 'Color', 'w');

    subplot(1, 2, 1);
    [x1, x2] = meshgrid(linspace(center(1) - 0.5, center(1) + 0.5, 100), ...
        linspace(center(2) - 0.5, center(2) + 0.5, 100));
    values = 0.5 * (H(1, 1) * x1.^2 + 2 * H(1, 2) * x1 .* x2 + ...
        H(2, 2) * x2.^2) + b(1) * x1 + b(2) * x2;
    contourf(x1, x2, values, 20, 'LineColor', 'none');
    hold on;
    sampleHandle = scatter(X(1, :), X(2, :), 18, 'k', 'filled');
    arrowLength = 0.28;
    trueHandle = quiver(center(1), center(2), ...
        arrowLength * trueDirection(1), arrowLength * trueDirection(2), ...
        0, 'Color', [0.1 0.65 0.2], 'LineWidth', 2, 'MaxHeadSize', 0.8);
    rawDirection = representative.noisyG / norm(representative.noisyG);
    rawHandle = quiver(center(1), center(2), ...
        arrowLength * rawDirection(1), arrowLength * rawDirection(2), ...
        0, 'Color', [0.85 0.2 0.15], 'LineWidth', 2, 'MaxHeadSize', 0.8);
    denoisedHandle = quiver(center(1), center(2), ...
        arrowLength * representative.denoised(1), ...
        arrowLength * representative.denoised(2), ...
        0, 'Color', [0.1 0.35 0.9], 'LineWidth', 2, 'MaxHeadSize', 0.8);
    hold off;
    axis equal;
    grid on;
    xlabel('x_1');
    ylabel('x_2');
    title('One noisy sample cloud and its directions');
    legend([sampleHandle, trueHandle, rawHandle, denoisedHandle], ...
        {'Sample locations', 'True direction', 'Noisy gradient', ...
         'Denoised direction'}, 'Location', 'best');
    colorbar;

    subplot(1, 2, 2);
    rawSorted = sort(rawErrors);
    denoisedSorted = sort(denoisedErrors);
    stairs(rawSorted, (1:numel(rawSorted))' / numel(rawSorted), ...
        'Color', [0.85 0.2 0.15], 'LineWidth', 2);
    hold on;
    stairs(denoisedSorted, (1:numel(denoisedSorted))' / numel(denoisedSorted), ...
        'Color', [0.1 0.35 0.9], 'LineWidth', 2);
    hold off;
    grid on;
    xlim([0 180]);
    ylim([0 1]);
    xlabel('Angular error from true direction (degrees)');
    ylabel('Fraction of trials at or below this error');
    title('Repeated trials: lower errors are better');
    legend({'Noisy gradient', 'Denoised direction'}, 'Location', 'southeast');
end
