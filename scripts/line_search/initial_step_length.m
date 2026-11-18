function alpha = initial_step_length(g, d, alphaMax)
% initial_step_length - Första steget i backtracking (line_search.m).
%
% Väljer alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF: brant lutning ger kort steg, flack lutning långt steg.
% Resultatet begränsas till [alphaMin, alphaMax]. Negativ lutning ger ett
% negativt värde och lutning 0 ger Inf; gränserna tar hand om båda.

    deltaF   = 1E2;   % förutsagd ökning i första punkten
    alphaMin = 0.25;  % kortaste första steg

    alpha = deltaF/(g'*d);
    alpha = min(alpha, alphaMax);
    alpha = max(alpha, alphaMin);
end
