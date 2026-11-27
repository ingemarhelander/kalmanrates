function d = find_sd_direction(x, g, xLbd, normalizeDirection)
% FIND_SD_DIRECTION Return the steepest-ascent direction, respecting lower bounds.

    if nargin < 4 || isempty(normalizeDirection)
        normalizeDirection = false;
    end

    d = g;
    if normalizeDirection && norm(d) > 0
        d = d / norm(d);
    end

    if g' * d < 0
        fprintf('This is not an ascent direction\n');
    end

    alphaLbd = ones(size(x)) * inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);

    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d * alphaLbdMin * 0.99;
    end

    if any(isnan(d))
        error('Search direction contains NaN');
    end
end
