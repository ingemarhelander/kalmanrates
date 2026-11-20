function alpha = initial_step_length_main(g, d, alphaMax, alphaMin)
% initial_step_length_main - Första steget i line_search_main.
%
% Väljer alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF, begränsat till [alphaMin, alphaMax] (alphaMin vinner om
% alphaMax < alphaMin; line_search_main klipper sedan mot alphaMax).

    deltaF = 100;     % förutsagd ökning i första punkten

    alpha = alphaMax;
    alphaObjInc = deltaF/(g'*d);
    
    alpha = min(alpha, alphaObjInc);
    alpha = max(alpha, alphaMin);

end
