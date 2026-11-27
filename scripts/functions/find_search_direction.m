function [d, D, needsDenoising] = find_search_direction(x, g, xLbd, D)
% find_search-direction Computes an orthonormalised search direction 
% using a modified conjugate gradient / Gram-Schmidt approach.
needsDenoising = false;

    %Initial direction
    d = g;
    
    if isempty(D)
        % First iteration: use normed gradient as first serach direction
        norm_d = norm(d);
        if norm_d > 0
            d = d / norm_d;
        end
    else
        % 2. Orthonormalisation
        % Gram-Schmidt
        U = orth(D);
        d = d - U * (U' * d);
        
        % Normalising new search direction
        norm_d = norm(d);
        if norm_d > 1e-10 * norm(g)
            d = d / norm_d;
        else
            if norm(g) > 0
               d = g / norm(g);
               needsDenoising = true;
            end
        end
    end

    % --- Bound checks ---
    
    % Largest step length (along d) before each variable hits its lower bound
    alphaLbd = ones(size(x))*inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);

    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d * alphaLbdMin * 0.99;   % Ensures that alpha = 1 is feasible
    end
    
    if any(isnan(d))
        error('Search direction contains NaN');
    end
end