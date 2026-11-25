function [extendedFit, extendedEvals] = extend_polynomial_fit( ...
    x, d, K, mu, varInd, evals, opts)
% Experimentell utvidgning av det senaste ML-polynomfitet.
% Ändrar inte optimeringens x, f, g eller valda steg.

    if nargin < 7 || isempty(opts)
        opts = struct();
    end
    opts = default_extension_opts(opts);

    if isempty(evals.fit)
        error('extend_polynomial_fit:NoFit', ...
              'evals.fit är tomt. Det finns inget polynom att utvidga.');
    end
    if exist('valid_line_search_point', 'file') ~= 2
        error('extend_polynomial_fit:MissingValidator', ...
              'valid_line_search_point.m kunde inte hittas på MATLAB-sökvägen.');
    end
    if ~isfield(evals, 'A') || ~isfield(evals, 'G')
        error('extend_polynomial_fit:MissingData', ...
              'evals måste innehålla fälten A och G.');
    end

    extendedEvals = evals;
    extendedFit = extendedEvals.fit{end};
    alphaMax = extendedEvals.fitAlphaMax{end};

    if ~(isscalar(alphaMax) && isfinite(alphaMax) && alphaMax > 0)
        error('extend_polynomial_fit:InvalidAlphaMax', ...
              'alphaMax måste vara ett ändligt positivt tal.');
    end

    alphaCurrent = 0;
    alphaStart = opts.startAlpha;

    % Återanvänd startpunkten om den redan finns.
    existing = find(abs(extendedEvals.A - alphaStart) <= ...
                    opts.duplicateTol, 1);

    if ~isempty(existing)
        alphaCurrent = extendedEvals.A(existing);
        fprintf('  Extension: alpha = %.6f finns redan.\n', alphaCurrent);
    else
        [alphaTry, isValid] = find_valid_alpha( ...
            x, d, K, varInd, alphaCurrent, alphaStart, opts);

        if ~isValid
            extendedEvals.extensionReason = ...
                'Ingen giltig startpunkt kunde hittas vid alpha = 0.95.';
            fprintf('  Extension stoppad: %s\n', ...
                    extendedEvals.extensionReason);
            return
        end

        [fTry, gTry] = kalmanObj(x + alphaTry*d, K, mu, varInd);

        if ~(isfinite(fTry) && all(isfinite(gTry)))
            extendedEvals.extensionReason = ...
                'Kalman returnerade NaN/Inf vid startpunkten.';
            fprintf('  Extension stoppad: %s\n', ...
                    extendedEvals.extensionReason);
            return
        end

        extendedEvals = append_evaluation( ...
            extendedEvals, alphaTry, x + alphaTry*d, fTry, gTry);

        alphaCurrent = alphaTry;

        [extendedFit, extendedEvals] = refit_polynomial( ...
            extendedEvals, alphaMax, d);

        fprintf('  Extension: startpunkt alpha = %.6f\n', alphaCurrent);
    end

    nNewEvals = 0;
    stopReason = '';

    while nNewEvals < opts.maxEvals

        tCurrent = alphaCurrent / alphaMax;
        derivative = poly_slope(extendedFit, tCurrent) / alphaMax;

        fprintf('  Extension: alpha = %.6f, dF/dalpha = %+.6e\n', ...
                alphaCurrent, derivative);

        if abs(derivative) >= opts.derivativeThreshold
            stopReason = sprintf( ...
                '|dF/dalpha| nådde %.3e.', opts.derivativeThreshold);
            break
        end

        deltaAlpha = opts.stepScale / ...
            max(abs(derivative), opts.slopeFloor);
        deltaAlpha = min(max(deltaAlpha, opts.minStep), opts.maxStep);

        alphaTarget = min(alphaCurrent + deltaAlpha, opts.alphaLimit);

        if alphaTarget <= alphaCurrent + opts.duplicateTol
            stopReason = 'Ingen ytterligare steglängd återstår före alphaLimit.';
            break
        end

        [alphaTry, isValid] = find_valid_alpha( ...
            x, d, K, varInd, alphaCurrent, alphaTarget, opts);

        if ~isValid
            stopReason = 'Ingen giltig punkt hittades över minStep.';
            break
        end

        existing = find(abs(extendedEvals.A - alphaTry) <= ...
                        opts.duplicateTol, 1);
        if ~isempty(existing)
            alphaCurrent = extendedEvals.A(existing);
            continue
        end

        [fTry, gTry] = kalmanObj(x + alphaTry*d, K, mu, varInd);

        % Backtracking även om Kalman trots giltighetskontrollen ger NaN/Inf.
        if ~(isfinite(fTry) && all(isfinite(gTry)))
            delta = alphaTry - alphaCurrent;
            accepted = false;

            while delta/2 >= opts.minStep
                delta = delta/2;
                alphaBack = alphaCurrent + delta;

                if ~valid_line_search_point(x + alphaBack*d, K, varInd)
                    continue
                end

                [fBack, gBack] = kalmanObj( ...
                    x + alphaBack*d, K, mu, varInd);

                if isfinite(fBack) && all(isfinite(gBack))
                    alphaTry = alphaBack;
                    fTry = fBack;
                    gTry = gBack;
                    accepted = true;
                    break
                end
            end

            if ~accepted
                stopReason = 'Kalman gav NaN/Inf och backtracking misslyckades.';
                break
            end
        end

        extendedEvals = append_evaluation( ...
            extendedEvals, alphaTry, x + alphaTry*d, fTry, gTry);

        alphaCurrent = alphaTry;
        nNewEvals = nNewEvals + 1;

        [extendedFit, extendedEvals] = refit_polynomial( ...
            extendedEvals, alphaMax, d);

        fprintf('  Extension: ny punkt alpha = %.6f, F = %.8f\n', ...
                alphaCurrent, fTry);
    end

    if isempty(stopReason)
        stopReason = sprintf('maxEvals = %d nåddes.', opts.maxEvals);
    end

    extendedEvals.extensionReason = stopReason;
    extendedEvals.extensionEvaluations = nNewEvals;

    fprintf('  Extension avslutad efter %d nya utvärderingar: %s\n', ...
            nNewEvals, stopReason);
end


function [alphaTry, isValid] = find_valid_alpha( ...
    x, d, K, varInd, alphaCurrent, alphaTarget, opts)

    delta = alphaTarget - alphaCurrent;
    isValid = false;
    alphaTry = NaN;

    while delta >= opts.minStep
        candidate = alphaCurrent + delta;

        if valid_line_search_point(x + candidate*d, K, varInd)
            alphaTry = candidate;
            isValid = true;
            return
        end

        delta = delta/2;
    end
end


function evals = append_evaluation(evals, alpha, xNew, fNew, gNew)

    evals.X(:, end+1) = xNew;
    evals.F(end+1) = fNew;
    evals.G(:, end+1) = gNew;
    evals.A(end+1) = alpha;
    evals.Q(end+1) = NaN;
    evals.AHat(end+1) = NaN;
end


function [fit, evals] = refit_polynomial(evals, alphaMax, d)

    t = evals.A(:) / alphaMax;
    F = evals.F(:);
    S = alphaMax * (evals.G' * d);

    [q, fit, info] = select_poly_degree(t, F, S);

    evals.fit{end+1} = fit;
    evals.fitInfo{end+1} = info;
    evals.fitAlpha{end+1} = t;
    evals.fitF{end+1} = F;
    evals.fitSlope{end+1} = S;
    evals.fitAlphaMax{end+1} = alphaMax;

    evals.Q(end) = q;
end


function opts = default_extension_opts(opts)

    def = struct( ...
        'startAlpha', 0.95, ...
        'alphaLimit', 1.3, ...
        'derivativeThreshold', 1e6, ...
        'stepScale', 0.02, ...
        'slopeFloor', 1e-8, ...
        'minStep', 1e-6, ...
        'maxStep', 0.05, ...
        'maxEvals', 50, ...
        'duplicateTol', 1e-10);

    unknown = setdiff(fieldnames(opts), fieldnames(def));
    if ~isempty(unknown)
        error('extend_polynomial_fit:UnknownOption', ...
              'Okänd inställning: %s', strjoin(unknown, ', '));
    end

    names = fieldnames(def);
    for k = 1:numel(names)
        if ~isfield(opts, names{k})
            opts.(names{k}) = def.(names{k});
        end
    end
end
