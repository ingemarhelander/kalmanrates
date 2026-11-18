function alpha = initial_step_length(g, d, alphaMax, alphaMin)
% initial_step_length - Första steget i linjesökningen (line_search_main
% och backtracking i line_search.m).
%
% Väljer alpha så att tangenten f + alpha*(g'*d) förutsäger ökningen
% deltaF: brant lutning ger kort steg, flack lutning långt steg.
% Resultatet begränsas till [alphaMin, alphaMax]; är alphaMax < alphaMin
% gäller alphaMax. Negativ lutning ger ett negativt värde och lutning 0 ger
% Inf; gränserna tar hand om båda.
%
%   alphaMax  längsta första steg   (standard 1)
%   alphaMin  kortaste första steg  (standard 1E-5)

    if nargin < 3 || isempty(alphaMax)
        alphaMax = 1;
    end
    if nargin < 4 || isempty(alphaMin)
        alphaMin = 1E-5;
    end

    deltaF = 100;     % förutsagd ökning i första punkten

    alpha = deltaF/(g'*d);
    alpha = max(alpha, alphaMin);
    alpha = min(alpha, alphaMax);
end
