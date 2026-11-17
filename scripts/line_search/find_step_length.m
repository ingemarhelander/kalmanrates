function alphaNew = find_step_length(fit, h, alphaMax)
% find_step_length - Nästa steglängd: polynomets första maximum längs linjen.
%
% Går från alpha = 0 längs ML-polynomet med fasta steg h (h, 2h, 3h, ...
% upp till alphaMax) och stannar första gången derivatan p'(alpha) blir
% negativ. Då har polynomet passerat ett lokalt maximum mellan de två sista
% stegen, och mittpunkten mellan dem returneras. Noggrannheten är h/2.
%
% Bara polynomet används, ingen kalmanObj, så h kan vara mycket litet.
% Blir p' aldrig negativ före alphaMax ligger polynomets maximum på
% gränsen, och alphaMax returneras.
%
% In:
%   fit       anpassningen från ml_poly_main
%   h         steglängd för promenaden längs polynomet
%   alphaMax  längsta tillåtna steg
% Ut:
%   alphaNew  nästa steglängd att utvärdera med kalmanObj

    % Alla steg på en gång (samma resultat som att gå ett steg i taget)
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
