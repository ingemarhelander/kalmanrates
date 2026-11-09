function d = find_search_direction(x, g, xLbd, D)
% FIND_SEARCH_DIRECTION Computes an orthonormalised search direction 
% using a modified conjugate gradient / Gram-Schmidt approach.

    % 1. Initial direction[cite: 1]
    d = g;
    
    if isempty(D)
        % Första iterationen: Gradienten väljs som initial sökriktning[cite: 1].
        norm_d = norm(d);
        if norm_d > 0
            d = d / norm_d;
        end
        fprintf('[find_search_direction] Använder normaliserad gradient som initial sökriktning.\n');
    else
        % 2. Orthonormalisation[cite: 1]
        % Gram-Schmidt ortogonalisering mot alla tidigare riktningar[cite: 1].
        fprintf('[find_search_direction] Utför Gram-Schmidt mot %d tidigare riktning(ar)...\n', size(D, 2));
        for i = 1:size(D, 2)
            d_prev = D(:, i);
            % Projicera d på d_prev och subtrahera
            proj = (d' * d_prev) * d_prev;
            d = d - proj;
        end
        
        % Normalisera den nya riktningen för att säkerställa att den är ortonormal[cite: 1]
        norm_d = norm(d);
        if norm_d > 1e-12
            d = d / norm_d;
            fprintf('[find_search_direction] Riktning framgångsrikt ortogonaliserad och normaliserad.\n');
        else
            fprintf('[find_search_direction] VARNING: Gradienten ryms helt i tidigare sökriktningar.\n');
            d = g; 
        end
    end

    % --- Originalkod för begränsningskontroll (Bounds) ---
    
    % Kontrollera ascent direction
    ascent_check = g' * d;
    if ascent_check < 0
        fprintf('[find_search_direction] VARNING: Detta är inte en ascent direction!\n');
    end

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