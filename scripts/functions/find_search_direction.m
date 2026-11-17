function [d, D] = find_search_direction(x, g, xLbd, D)
% find_search-direction Computes an orthonormalised search direction 
% using a modified conjugate gradient / Gram-Schmidt approach.

    % Aktiv mängd: variabler som ligger vid sin undre gräns och som
    % gradienten vill minska låses (deras komponent sätts till 0). Annars
    % blir alphaLbdMin ~ 0 nedan och hela d skalas ner till noll.
    atBound = (xLbd ~= -inf) & (x - xLbd <= 1e-8*max(1, abs(xLbd))) & (g < 0);
    g(atBound) = 0;
    if any(atBound)
        fprintf('[find_search_direction] %d variabel/variabler vid undre gränsen låses.\n', nnz(atBound));
    end

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
        % Gram-Schmidt ortogonalisering mot alla tidigare riktningar
        fprintf('[find_search_direction] Utför Gram-Schmidt mot %d tidigare riktning(ar)...\n', size(D, 2));
        % Kolumnerna i D är skalade mot gränserna och inte ortonormala,
        % så projicera mot en ortonormal bas för span(D)
        U = orth(D);
        d = d - U * (U' * d);
        d(atBound) = 0;   % projektionen kan ha lagt tillbaka komponenter där

        % Ortonoaliseringkontroll
        orth_error = max(abs(U' * d));
        fprintf('[find_search_direction] [KONTROLL] Max ortogonaliseringsfel (U''*d): %.3e\n', orth_error);
        if orth_error > 1e-6
            fprintf('[find_search_direction] VARNING: Dålig ortogonalisering upptäcktes!\n');
        end

        % Normalisera den nya riktningen för att säkerställa att den är ortonormal
        norm_d = norm(d);
        if norm_d > 1e-10 * norm(g)
            d = d / norm_d;
            % Normeringskorntroll
            final_norm = norm(d);
            fprintf('[find_search_direction] [KONTROLL] Slutlig norm(d) efter normalisering: %.6f (bör vara 1.0)\n', final_norm);
            if abs(final_norm - 1.0) > 1e-6
                fprintf('[find_search_direction] VARNING: Normaliseringen misslyckades!\n');
            end

            fprintf('[find_search_direction] Riktning framgångsrikt ortogonaliserad och normaliserad.\n');
        else
            fprintf('[find_search_direction] VARNING: Gradienten ryms helt i tidigare sökriktningar.\n');
            if norm(g) > 0
               d = g / norm(g);
               D = [];
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
    % Diagnos: vilken variabel begränsar steget (kan tas bort senare)
    [m, i] = min(alphaLbd);
    fprintf('[find_search_direction] alphaLbdMin = %.3e (var %d, x-xLbd = %.3e, d = %.3e)\n', ...
            m, i, x(i) - xLbd(i), d(i));
    if alphaLbdMin < 1
        d = d * alphaLbdMin * 0.99;   % Ensures that alpha = 1 is feasible
    end

    if any(isnan(d))
        error('Search direction contains NaN');
    end
end