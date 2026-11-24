function alpha = initial_step_length_main(g, d, alphaMin)
% initial_step_length_main - Första steget i line_search_main.
%
% Väljer alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF, begränsat till [alphaMin, 1]. alpha = 1 är alltid tillåtet,
% eftersom riktningen d redan är skalad mot variansgränserna.

    deltaF = 100;     % förutsagd ökning i första punkten

    alpha = 1;
    alphaObjInc = deltaF/(g'*d);

    alpha = min(alpha, alphaObjInc);
    alpha = max(alpha, alphaMin);

end
