function stop = evaluate(fit, alphaHat, alphas, ok)
% evaluate - Ska line search stanna innan alphaHat utvärderas?
%
% TILLFÄLLIG VERSION: innehåller bara stoppvillkoren som tidigare låg i
% line_search_main. Villkoren i diagrammet är inte implementerade än:
%   - ny gradient mindre än ½ gammal?
%   - ger alpha en betydande förbättring?
%   - p(alpha) ...
%
% In:
%   fit       ML-anpassningen från ml_poly_main
%   alphaHat  föreslagen steglängd från find_step_length
%   alphas    steglängder för alla utvärderade punkter (startpunkten 0 först)
%   ok        true för punkter med ändliga F och G
% Ut:
%   stop      true om line search ska stanna

    dMinRel = 1e-3;   % minsta avstånd till en utvärderad punkt, relativt största steget (Req. 31)

    % För nära en utvärderad punkt (Req. 31). Relativt stegen och inte
    % alphaMax, som kan vara mycket större än stegen.
    tooClose = min(abs(alphaHat - alphas)) < dMinRel*max(alphas);

    % Förutsagd ökning mindre än bruset. fit.sigmaF är ~0 vid exakt
    % anpassning, då stoppar bara de andra kriterierna.
    noGain = poly_value(fit, alphaHat) - max(poly_value(fit, alphas(ok))) < fit.sigmaF;

    stop = tooClose || noGain;
end
