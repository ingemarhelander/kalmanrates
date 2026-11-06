function d = find_search_direction(x, g, xLbd, pastDirections)
% Steepest ascent direction, shortened so that a full step (alpha = 1)
% stays inside the lower bounds and is orthogonal to pastDirections.

    if nargin < 4 || isempty(pastDirections)
        d = g;
    else
        if size(pastDirections, 1) ~= numel(g)
            error('find_search_direction:DimensionMismatch', ...
                'Past directions must have the same number of rows as g has elements.');
        end

        % Build an orthonormal basis for the span of the past directions.
        basis = zeros(size(pastDirections));
        basisSize = 0;
        for j = 1:size(pastDirections, 2)
            v = pastDirections(:, j);
            originalNorm = norm(v);
            if originalNorm == 0
                continue;
            end

            for pass = 1:2
                for k = 1:basisSize
                    v = v - basis(:, k) * (basis(:, k)' * v);
                end
            end

            residualNorm = norm(v);
            if residualNorm > sqrt(eps) * originalNorm
                basisSize = basisSize + 1;
                basis(:, basisSize) = v / residualNorm;
            end
        end
        basis = basis(:, 1:basisSize);

        % Remove the gradient component in the span of past directions.
        d = g;
        for pass = 1:2
            for j = 1:basisSize
                d = d - basis(:, j) * (basis(:, j)' * d);
            end
        end
    end

    if norm(d) <= sqrt(eps) * norm(g)
        d = zeros(size(g));
        return;
    end

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