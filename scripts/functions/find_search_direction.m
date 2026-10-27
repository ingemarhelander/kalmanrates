function d = find_search_direction(x, g, xLbd)
% Steepest ascent direction, shortened so that a full step (alpha = 1)
% stays inside the lower bounds.

    d = g;

    if g'*d < 0
        fprintf('This is not an ascent direction\n');
    end

    % Largest step length (along d) before each variable hits its lower bound
    alphaLbd = ones(size(x))*inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);

    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d*alphaLbdMin*0.99;   % Ensures that alpha = 1 is feasible
    end

    if any(isnan(d))
        error('Search direction contains NaN');
    end
end