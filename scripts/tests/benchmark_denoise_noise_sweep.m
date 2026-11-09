function [summary, results, figureHandle] = benchmark_denoise_noise_sweep()
%BENCHMARK_DENOISE_NOISE_SWEEP Quantify denoise performance under noise.
%   [SUMMARY, RESULTS, FIGURE] = benchmark_denoise_noise_sweep() runs a
%   reproducible Monte Carlo experiment on a ten-dimensional quadratic.
%   It sweeps sample counts around 40 and the measurement-noise level, prints
%   direction-error metrics, and opens a four-panel summary figure.
%
%   Gradient noise is reported as the RMS noise-vector norm divided by
%   norm(trueGradientAtCenter).
%   Function-value noise is scaled to the local objective variation:
%   sigma_f = noiseLevel * norm(trueGradientAtCenter) * cloudRadius.
%   The reported "noise removed" is the reduction in direction RMSE relative
%   to using the noisy gradient at the center; denoise returns a unit search
%   direction, so this does not claim to recover the measurement noise's
%   absolute magnitude.

    addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'functions'));

    sampleCounts = [25, 32, 36, 40, 44, 48, 55];
    noiseLevels = [0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8];
    trialCount = 20;
    cloudRadius = 0.12;
    successAngleDeg = 10;
    requiredSuccessRate = 0.80;
    rng(20261109, 'twister');

    dimension = 10;
    randomMatrix = randn(dimension);
    H = eye(dimension) + (randomMatrix' * randomMatrix) / dimension;
    center = linspace(-0.3, 0.3, dimension)';
    b = linspace(0.6, -0.4, dimension)';
    trueGradient = H * center + b;
    trueDirection = trueGradient / norm(trueGradient);

    numberOfCounts = numel(sampleCounts);
    numberOfNoiseLevels = numel(noiseLevels);
    directionErrorRaw = nan(numberOfCounts, numberOfNoiseLevels, trialCount);
    directionErrorDenoised = nan(numberOfCounts, numberOfNoiseLevels, trialCount);
    fitSucceeded = false(numberOfCounts, numberOfNoiseLevels, trialCount);

    for countIndex = 1:numberOfCounts
        sampleCount = sampleCounts(countIndex);
        for noiseIndex = 1:numberOfNoiseLevels
            noiseLevel = noiseLevels(noiseIndex);
            % Per-component standard deviation gives the requested RMS
            % gradient-noise norm after combining all dimensions.
            sigmaGradient = noiseLevel * norm(trueGradient) / sqrt(dimension);
            sigmaFunction = noiseLevel * norm(trueGradient) * cloudRadius;

            for trial = 1:trialCount
                X = [center + cloudRadius * randn(dimension, sampleCount - 1), center];
                F = zeros(1, sampleCount);
                G = zeros(dimension, sampleCount);

                for sample = 1:sampleCount
                    point = X(:, sample);
                    F(sample) = 0.5 * point' * H * point + b' * point;
                    G(:, sample) = H * point + b;
                end

                noisyF = F + sigmaFunction * randn(size(F));
                noisyG = G + sigmaGradient * randn(size(G));

                rawDirection = noisyG(:, end) / max(norm(noisyG(:, end)), eps);
                directionErrorRaw(countIndex, noiseIndex, trial) = ...
                    directionAngleDeg(rawDirection, trueDirection);

                evals.X = X;
                evals.F = noisyF;
                evals.G = noisyG;

                try
                    % denoise prints diagnostic messages on every call; capture
                    % them so the benchmark reports only its useful summaries.
                    evalc(['denoisedDirection = denoise(center, noisyF(end), ', ...
                        'noisyG(:, end), eye(dimension), {evals}, sampleCount, ', ...
                        '[], [], [], -inf(dimension, 1));']);
                    if all(isfinite(denoisedDirection)) && norm(denoisedDirection) > 0
                        directionErrorDenoised(countIndex, noiseIndex, trial) = ...
                            directionAngleDeg(denoisedDirection, trueDirection);
                        fitSucceeded(countIndex, noiseIndex, trial) = true;
                    end
                catch
                    % Keep failed fits in the success-rate denominator.
                end
            end
        end
    end

    rawRmseDeg = nan(numberOfCounts, numberOfNoiseLevels);
    denoisedRmseDeg = nan(numberOfCounts, numberOfNoiseLevels);
    reductionPercent = nan(numberOfCounts, numberOfNoiseLevels);
    successRatePercent = zeros(numberOfCounts, numberOfNoiseLevels);
    angleSuccessRatePercent = zeros(numberOfCounts, numberOfNoiseLevels);
    failureRatePercent = zeros(numberOfCounts, numberOfNoiseLevels);

    for countIndex = 1:numberOfCounts
        for noiseIndex = 1:numberOfNoiseLevels
            rawErrors = reshape(directionErrorRaw(countIndex, noiseIndex, :), 1, []);
            denoisedErrors = reshape(directionErrorDenoised(countIndex, noiseIndex, :), 1, []);
            rawRmseDeg(countIndex, noiseIndex) = sqrt(mean(rawErrors .^ 2, 'omitnan'));
            denoisedRmseDeg(countIndex, noiseIndex) = ...
                sqrt(mean(denoisedErrors .^ 2, 'omitnan'));
            successRatePercent(countIndex, noiseIndex) = ...
                100 * mean(reshape(fitSucceeded(countIndex, noiseIndex, :), 1, []));
            angleSuccessRatePercent(countIndex, noiseIndex) = 100 * mean( ...
                reshape(fitSucceeded(countIndex, noiseIndex, :), 1, []) ...
                & denoisedErrors <= successAngleDeg);
            failureRatePercent(countIndex, noiseIndex) = 100 - ...
                successRatePercent(countIndex, noiseIndex);

            if rawRmseDeg(countIndex, noiseIndex) > 100 * eps
                reductionPercent(countIndex, noiseIndex) = 100 * ...
                    (1 - denoisedRmseDeg(countIndex, noiseIndex) / ...
                    rawRmseDeg(countIndex, noiseIndex));
            end
        end
    end

    [sampleGrid, noiseGrid] = ndgrid(sampleCounts, noiseLevels);
    summary = table(sampleGrid(:), noiseGrid(:), rawRmseDeg(:), ...
        denoisedRmseDeg(:), reductionPercent(:), successRatePercent(:), ...
        angleSuccessRatePercent(:), ...
        failureRatePercent(:), 'VariableNames', ...
        {'SamplePoints', 'RMSGradientNoiseOverGradientNorm', 'RawDirectionRMSEdeg', ...
        'DenoisedDirectionRMSEdeg', 'DirectionRMERReductionPercent', ...
        'FitSuccessPercent', 'AngleSuccessPercent', 'FitFailurePercent'});

    fprintf('\nDenoise benchmark: 10D quadratic model, %d Monte Carlo trials/configuration.\n', ...
        trialCount);
    fprintf('Noise level = RMS gradient-noise norm / ||g(center)||; success = direction error <= %.1f deg.\n', ...
        successAngleDeg);
    disp(summary);

    requiredCounts = nan(1, numberOfNoiseLevels);
    for noiseIndex = 1:numberOfNoiseLevels
        successful = angleSuccessRatePercent(:, noiseIndex) >= 100 * requiredSuccessRate ...
            & denoisedRmseDeg(:, noiseIndex) <= successAngleDeg ...
            & reductionPercent(:, noiseIndex) > 0;
        firstSuccessful = find(successful, 1, 'first');
        if ~isempty(firstSuccessful)
            requiredCounts(noiseIndex) = sampleCounts(firstSuccessful);
        end
    end

    fprintf('\nMinimum tested sample count meeting all criteria (>=%.0f%% trials below angle target, ', ...
        100 * requiredSuccessRate);
    fprintf('RMSE <= %.1f deg, and positive RMSE reduction):\n', successAngleDeg);
    for noiseIndex = 1:numberOfNoiseLevels
        if noiseLevels(noiseIndex) == 0
            fprintf('  noise %.2f: not applicable (no injected noise)\n', ...
                noiseLevels(noiseIndex));
        elseif isnan(requiredCounts(noiseIndex))
            fprintf('  noise %.2f: not reached in tested range\n', noiseLevels(noiseIndex));
        else
            fprintf('  noise %.2f: %d points\n', noiseLevels(noiseIndex), ...
                requiredCounts(noiseIndex));
        end
    end

    figureHandle = figure('Name', 'Denoise noise and sample-count benchmark', ...
        'Color', 'w', 'Position', [100, 100, 1200, 820]);
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, sprintf('10D quadratic denoising benchmark (%d trials per setting)', ...
        trialCount), 'FontSize', 15, 'FontWeight', 'bold');

    % Panel 1: sensitivity to noise at representative sample counts.
    ax1 = nexttile(layout);
    hold(ax1, 'on');
    representativeCountIndices = unique([1, 3, 5, numberOfCounts]);
    for countIndex = representativeCountIndices
        plot(ax1, noiseLevels, denoisedRmseDeg(countIndex, :), '-o', ...
            'LineWidth', 1.7, 'MarkerSize', 5, ...
            'DisplayName', sprintf('%d points', sampleCounts(countIndex)));
    end
    plot(ax1, noiseLevels, rawRmseDeg(end, :), 'k--', 'LineWidth', 1.6, ...
        'DisplayName', 'Raw center gradient');
    yline(ax1, successAngleDeg, ':', '10^\circ target', 'HandleVisibility', 'off');
    grid(ax1, 'on');
    xlabel(ax1, 'RMS gradient noise / ||g(x_0)||');
    ylabel(ax1, 'Direction RMSE (degrees)');
    title(ax1, 'Noise tolerance');
    legend(ax1, 'Location', 'northwest');

    % Panel 2: convergence as more sample points are available.
    ax2 = nexttile(layout);
    hold(ax2, 'on');
    selectedNoise = [0.1, 0.4, 0.8];
    for noiseLevel = selectedNoise
        noiseIndex = find(abs(noiseLevels - noiseLevel) < 1e-12, 1);
        plot(ax2, sampleCounts, denoisedRmseDeg(:, noiseIndex), '-o', ...
            'LineWidth', 1.7, 'MarkerSize', 5, ...
            'DisplayName', sprintf('Noise %.1f', noiseLevel));
    end
    yline(ax2, successAngleDeg, ':', '10^\circ target', 'HandleVisibility', 'off');
    grid(ax2, 'on');
    xlabel(ax2, 'Number of data points (including center)');
    ylabel(ax2, 'Denoised direction RMSE (degrees)');
    title(ax2, 'Effect of sample count');
    legend(ax2, 'Location', 'northeast');

    % Panel 3: chance of meeting the angle threshold, including fit failures.
    ax3 = nexttile(layout);
    imagesc(ax3, sampleCounts, noiseLevels, angleSuccessRatePercent');
    set(ax3, 'YDir', 'normal');
    colormap(ax3, parula(100));
    colorbar(ax3);
    caxis(ax3, [0, 100]);
    hold(ax3, 'on');
    contour(ax3, sampleCounts, noiseLevels, angleSuccessRatePercent', ...
        [100 * requiredSuccessRate, 100 * requiredSuccessRate], ...
        'k--', 'LineWidth', 1.5);
    grid(ax3, 'on');
    xlabel(ax3, 'Number of data points');
    ylabel(ax3, 'RMS gradient noise / ||g(x_0)||');
    title(ax3, 'Trials with direction error <= 10^\circ');

    % Panel 4: direction-error reduction compared with the raw center gradient.
    ax4 = nexttile(layout);
    hold(ax4, 'on');
    for noiseLevel = selectedNoise
        noiseIndex = find(abs(noiseLevels - noiseLevel) < 1e-12, 1);
        plot(ax4, sampleCounts, reductionPercent(:, noiseIndex), '-o', ...
            'LineWidth', 1.7, 'MarkerSize', 5, ...
            'DisplayName', sprintf('Noise %.1f', noiseLevel));
    end
    yline(ax4, 0, 'k:', 'No improvement', 'HandleVisibility', 'off');
    grid(ax4, 'on');
    xlabel(ax4, 'Number of data points');
    ylabel(ax4, 'Direction RMSE reduction (%)');
    title(ax4, 'How much center-direction error is removed');
    legend(ax4, 'Location', 'best');

    results.sampleCounts = sampleCounts;
    results.dimension = dimension;
    results.noiseLevels = noiseLevels;
    results.trialCount = trialCount;
    results.successAngleDeg = successAngleDeg;
    results.requiredSuccessRate = requiredSuccessRate;
    results.rawDirectionErrorDeg = directionErrorRaw;
    results.denoisedDirectionErrorDeg = directionErrorDenoised;
    results.fitSucceeded = fitSucceeded;
    results.rawDirectionRMSEDeg = rawRmseDeg;
    results.denoisedDirectionRMSEDeg = denoisedRmseDeg;
    results.directionRMERReductionPercent = reductionPercent;
    results.fitSuccessPercent = successRatePercent;
    results.angleSuccessPercent = angleSuccessRatePercent;
    results.fitFailurePercent = failureRatePercent;
    results.minimumSamplePointsForCriteria = requiredCounts;
end

function angle = directionAngleDeg(direction, reference)
    direction = direction / max(norm(direction), eps);
    cosine = max(-1, min(1, dot(direction, reference)));
    angle = acosd(cosine);
end