function alpha = initial_step_length(g, d, alphaMax)
% initial_step_length - Första steget i linjesökningen (line_search och
% line_search_main).
%
% Väljer alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF: brant lutning ger kort steg, flack lutning långt steg.
% Resultatet är högst alphaMax (lutning 0 ger Inf, som då blir alphaMax).

    deltaF = 100;     % förutsagd ökning i första punkten
    alpha = min(deltaF/(g'*d), alphaMax);
end
