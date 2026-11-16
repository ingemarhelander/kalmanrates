function [q, fit, info] = select_poly_degree(fits, method, options)
% SELECT_POLY_DEGREE Select a polynomial degree from precomputed fits.
%
%   [q, fit, info] = select_poly_degree(fits)
%   [q, fit, info] = select_poly_degree(fits, method)
%   [q, fit, info] = select_poly_degree(fits, method, options)
%
% INPUTS
%   fits    Cell array (or struct array) returned by mult_q_poly_fits.
%           Each fit must contain:
%             sigmaF, sigmaG, rhoF, rhoG, pTrendF, pTrendG, converged
%   method  'threshold' (default) or 'score'.
%   options Optional struct overriding the defaults below.
%
% OUTPUTS
%   q       Selected polynomial degree. Empty if threshold method finds
%           no acceptable degree.
%   fit     The selected fit, or [] if no degree was selected.
%   info    Diagnostics, acceptance flags and (for score method) scores.
%
% IMPORTANT
%   This function does NOT fit any models. It only selects among the
%   models already calculated by mult_q_poly_fits.
%
% METHOD 'threshold'
%   Selects the lowest degree satisfying all diagnostic thresholds.
%   Right now: abs(rhoF) <=0.2
%              abs(rhoG) <=0.2
%              ptrendF >= 0.05
%              ptrendG >= 0.05
%              converged == true
%
% METHOD 'score'
%   Scores the fits using residual variance, residual autocorrelation,
%   trend-test p-values and model complexity. The lowest score wins.
%              1. Residual variance (sigmaF, sigmaG)
%              2. Residual autocorrelation (rhoF, rhoG)
%              3. Trend-test p-values (pTrendF, pTrendG)
%              4. Model complexity (degree)
%
% DEFAULT OPTIONS
%   options.weights         = [0.35 0.30 0.25 0.10]
%       Weights correspond to variance, autocorrelation, trend and
%       complexity, respectively. They are normalized to sum to one.

    if nargin < 2 || isempty(method)
        method = 'threshold';
    end
    if nargin < 3 || isempty(options)
        options = struct();
    end

    if isstring(method)
        method = char(method);
    end
    if ~ischar(method)
        error('method must be ''threshold'' or ''score''.');
    end
    method = lower(strtrim(method));

    if ~isstruct(options)
        error('options must be a struct.');
    end

    % Default settings (can be changed)
    defaults.maxAbsRho = 0.20;
    defaults.minTrendP = 0.05;
    defaults.weights = [0.35 0.30 0.25 0.10];

    names = fieldnames(defaults);
    for k = 1:numel(names)
        name = names{k};
        if ~isfield(options, name) || isempty(options.(name))
            options.(name) = defaults.(name);
        end
    end

    if ~isscalar(options.maxAbsRho) || options.maxAbsRho < 0
        error('options.maxAbsRho must be a nonnegative scalar.');
    end
    if ~isscalar(options.minTrendP) || options.minTrendP < 0 || ...
            options.minTrendP > 1
        error('options.minTrendP must be between 0 and 1.');
    end
    if numel(options.weights) ~= 4 || any(options.weights < 0) || ...
            sum(options.weights) <= 0
        error('options.weights must contain four nonnegative values, not all zero.');
    end
    weights = options.weights(:)' / sum(options.weights);

    % Convert struct array to cell array for uniform handling.
    if isstruct(fits)
        fits = num2cell(fits);
    end
    if ~iscell(fits) || isempty(fits)
        error('fits must be a nonempty cell array or struct array.');
    end

    nFits = numel(fits);
    required = {'sigmaF','sigmaG','rhoF','rhoG', ...
                'pTrendF','pTrendG','converged'};

    sigmaF = nan(nFits,1);
    sigmaG = nan(nFits,1);
    rhoF = nan(nFits,1);
    rhoG = nan(nFits,1);
    pTrendF = nan(nFits,1);
    pTrendG = nan(nFits,1);
    converged = false(nFits,1);

    for k = 1:nFits
        currentFit = fits{k};
        if ~isstruct(currentFit)
            error('fits{%d} must be a struct.', k);
        end

        for j = 1:numel(required)
            if ~isfield(currentFit, required{j})
                error('fits{%d} is missing the field "%s".', ...
                      k, required{j});
            end
        end

        sigmaF(k) = currentFit.sigmaF;
        sigmaG(k) = currentFit.sigmaG;
        rhoF(k) = currentFit.rhoF;
        rhoG(k) = currentFit.rhoG;
        pTrendF(k) = currentFit.pTrendF;
        pTrendG(k) = currentFit.pTrendG;
        converged(k) = logical(currentFit.converged);
    end

    degree = (1:nFits)';
    valid = converged & isfinite(sigmaF) & isfinite(sigmaG) & ...
            sigmaF > 0 & sigmaG > 0 & isfinite(rhoF) & ...
            isfinite(rhoG) & isfinite(pTrendF) & isfinite(pTrendG);

    % Method A: choose the lowest degree that passes every threshold.
    accepted = valid & abs(rhoF) <= options.maxAbsRho & ...
               abs(rhoG) <= options.maxAbsRho & ...
               pTrendF >= options.minTrendP & ...
               pTrendG >= options.minTrendP;

    % Initialize outputs.
    q = [];
    fit = [];
    scores = nan(nFits,1);

    switch method
        case 'threshold'
            selected = find(accepted, 1, 'first');
            if ~isempty(selected)
                q = degree(selected);
                fit = fits{selected};
                status = 'selected';
            else
                status = 'no_degree_passed';
            end

        case 'score'
            candidates = find(valid);
            if isempty(candidates)
                status = 'no_valid_fits';
            else
                % Lower values are better for all four score components.
                varianceRaw = 0.5 * ( ...
                    log(sigmaF(candidates) ./ min(sigmaF(candidates))) + ...
                    log(sigmaG(candidates) ./ min(sigmaG(candidates))) );

                autocorrRaw = 0.5 * (abs(rhoF(candidates)) + ...
                                     abs(rhoG(candidates)));

                % A p-value below minTrendP adds a penalty. No penalty
                % is added when the test does not indicate a trend.
                trendRaw = 0.5 * ( ...
                    max(0, log10(options.minTrendP ./ ...
                        max(pTrendF(candidates), realmin))) + ...
                    max(0, log10(options.minTrendP ./ ...
                        max(pTrendG(candidates), realmin))) );

                complexityRaw = degree(candidates) / max(degree(candidates));

                varianceScore = normalize_score(varianceRaw);
                autocorrScore = normalize_score(autocorrRaw);
                trendScore = normalize_score(trendRaw);
                complexityScore = normalize_score(complexityRaw);

                candidateScores = weights(1) * varianceScore + ...
                                  weights(2) * autocorrScore + ...
                                  weights(3) * trendScore + ...
                                  weights(4) * complexityScore;
                scores(candidates) = candidateScores;

                [~, localIndex] = min(candidateScores);
                selected = candidates(localIndex);
                q = degree(selected);
                fit = fits{selected};
                status = 'selected';
            end

        otherwise
            error('Unknown method "%s". Use ''threshold'' or ''score''.', method);
    end

    info.method = method;
    info.status = status;
    info.degree = degree;
    info.valid = valid;
    info.accepted = accepted;
    info.sigmaF = sigmaF;
    info.sigmaG = sigmaG;
    info.rhoF = rhoF;
    info.rhoG = rhoG;
    info.pTrendF = pTrendF;
    info.pTrendG = pTrendG;
    info.scores = scores;
    info.options = options;
    info.weights = weights;
end

function normalized = normalize_score(values)
% Map a score vector to [0,1]. A constant vector becomes all zeros.
    values = values(:);
    low = min(values);
    high = max(values);
    if high - low <= eps(max(1, max(abs(values))))
        normalized = zeros(size(values));
    else
        normalized = (values - low) / (high - low);
    end
end
