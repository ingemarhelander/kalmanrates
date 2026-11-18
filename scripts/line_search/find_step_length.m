function alphaNew = find_step_length(fit, h, alphaMax)
% find_step_length - Nästa steglängd: polynomets första maximum.
%
% Går längs polynomet från 0 med steg h (h, 2h, ... upp till alphaMax) och
% stannar första gången p'(alpha) < 0. Maximum ligger då mellan de två
% sista stegen, och mittpunkten returneras (noggrannhet h/2). Blir p'
% aldrig negativ returneras alphaMax. Ingen kalmanObj körs här.
%
% In:   fit       ML-anpassningen (från select_poly_degree)
%       h         steglängd längs polynomet
%       alphaMax  längsta tillåtna steg
% Ut:   alphaNew  nästa steglängd att köra kalmanObj i

    steps = min((1:ceil(alphaMax/h))*h, alphaMax);
    k = find(poly_slope(fit, steps) < 0, 1);

    if isempty(k)
        alphaNew = alphaMax;
        return
    end
    if k == 1
        prev = 0;
    else
        prev = steps(k - 1);
    end
    alphaNew = (prev + steps(k))/2;
end
