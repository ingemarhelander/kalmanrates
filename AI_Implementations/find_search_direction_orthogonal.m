function [d, D] = find_search_direction_orthogonal(x, g, xLbd, D)
% find_search_direction_orthogonal returns an ascent direction orthogonal to D.

    if isempty(D)
        d = g;
    else
        % Project into the null space of D' using an SVD-based decomposition.
        N = null(D');
        d = N * (N' * g);
    end

    norm_d = norm(d);
    if norm_d > 1e-10 * norm(g)
        d = d / norm_d;
    elseif norm(g) > 0
        % Match the original fallback when g lies in the span of D.
        d = g / norm(g);
        D = [];
    end

    % Scale the direction so a unit step does not cross a lower bound.
    alphaLbd = ones(size(x))*inf;
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
