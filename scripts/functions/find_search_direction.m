function [d, D] = find_search_direction(x, g, xLbd, D)
% find_search-direction Computes an orthonormalised search direction 
% using a modified conjugate gradient / Gram-Schmidt approach.

    %Initial direction
    d = g;
    messages = {};
    
    if isempty(D)
        % First iteration: use normed gradient as first serach direction
        norm_d = norm(d);
        if norm_d > 0
            d = d / norm_d;
        end
        %fprintf('[find_search_direction] Using normalized gradient as initial search direction.\n');
    else
        % 2. Orthonormalisation
        % Gram-Schmidt
        %fprintf('[find_search_direction]Performing Gram-Schmidt against %d previous direction(s)...\n', size(D, 2));
        U = orth(D);
        d = d - U * (U' * d);
        
        % Orthogonalisation check
        orth_error = max(abs(U' * d));
        %fprintf('[find_search_direction] [CHECK] Max orthogonalisation error (U''*d): %.3e\n', orth_error);
        if orth_error > 1e-6
            messages{end+1} = 'orthogonalization failed';
        end

        % Normalising new search direction
        norm_d = norm(d);
        if norm_d > 1e-10 * norm(g)
            d = d / norm_d;
           
            % Normaliseing check
            final_norm = norm(d);
            %fprintf('[find_search_direction] [CHECK] Final norm after normalisation: %.6f (should be 1.0)\n', final_norm);
            if abs(final_norm - 1.0) > 1e-6
                messages{end+1} = 'normalization failed';
            end

            %fprintf('[find_search_direction] Direction sucessfullly orthonormalised.\n');
        else
            messages{end+1} = 'gradient already lies in span(D)';
            if norm(g) > 0
               d = g / norm(g);
               D = [];
            end
        end
    end

    % --- Bound checks ---
    
    % Check ascent direction
    ascent_check = g' * d;
    if ascent_check <= 0
        messages{end+1} = sprintf( ...
            'direction is not ascent: g''*d=%+.3e, norm(g)=%.3e, norm(d)=%.3e', ...
            ascent_check, norm(g), norm(d));
    end

    % Largest step length (along d) before each variable hits its lower bound
    alphaLbd = ones(size(x))*inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);

    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d * alphaLbdMin * 0.99;   % Ensures that alpha = 1 is feasible
    end
    
    if ~isempty(messages)
        fprintf('[find_search_direction] WARNING: %s.\n', strjoin(messages, '; '));
    end

    if any(isnan(d))
        error('Search direction contains NaN');
    end
end