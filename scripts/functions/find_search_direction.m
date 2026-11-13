function d = find_search_direction(x, g, xLbd, D)
% find_search-direction Computes an orthonormalised search direction 
% using a modified conjugate gradient / Gram-Schmidt approach.

    %Initial direction
    d = g;
    
    if isempty(D)
        % Första iterationen: Gradienten väljs som initial sökriktning
        norm_d = norm(d);
        if norm_d > 0
            d = d / norm_d;
        end
        fprintf('[find_search_direction] Använder normaliserad gradient som initial sökriktning.\n');
    else
        % 2. Orthonormalisation
        % Ortogonalisera mot alla tidigare riktningar. Kolumnerna i D har
        % skalats mot de nedre gränserna (se nedan) och har därför inte
        % längd 1, och de behöver inte vara inbördes ortogonala. Därför
        % projiceras d bort från en ortonormal bas för span(D), och inte
        % mot kolumnerna en och en. Projektionen görs två gånger för att
        % ta bort avrundningsfel (klassisk Gram-Schmidt med omortogonalisering).
        fprintf('[find_search_direction] Utför Gram-Schmidt mot %d tidigare riktning(ar)...\n', size(D, 2));
        U = orth(D);
        for pass = 1:2
            d = d - U * (U' * d);
        end
        
        % Normalisera den nya riktningen för att säkerställa att den är ortonormal.
        % Toleransen är relativ till gradienten, som kan vara mycket stor.
        norm_d = norm(d);
        if norm_d > 1e-10 * norm(g)
            d = d / norm_d;
            fprintf('[find_search_direction] Riktning framgångsrikt ortogonaliserad och normaliserad (andel av gradienten kvar: %.2e).\n', ...
                    norm_d / norm(g));
        else
            fprintf('[find_search_direction] VARNING: Gradienten ryms helt i tidigare sökriktningar. Använder gradienten.\n');
            d = g;
            if norm(g) > 0
                d = g / norm(g);
            end
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