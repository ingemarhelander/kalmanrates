function alpha = initial_step_length_main(g, d, alphaMin, alphaMax, denoise)
% initial_step_length_main - Första steget i line_search_main.
%
% CG-steg: alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF. Denoise-steg: alphaMax (den uppmätta lutningen ska inte styra).
% Resultatet begränsas till [alphaMin, alphaMax].

    deltaF = 100;     % förutsagd ökning i första punkten

    alpha = alphaMax;
    alphaObjInc = deltaF/(g'*d);

    alpha = min(alpha, alphaObjInc);
    alpha = max(alpha, alphaMin);

end
