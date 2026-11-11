function [alphaNew, alphaMax] = find_step_length(fit, alphaPrev, alphaMax, isBetter)
% find_step_length - Nästa steglängd från ML-polynomet (Oskar).
%
% TILLFÄLLIG VERSION: söker polynomets max på ett rutnät, bara för att
% line_search_main ska gå att köra. Antar att fit.y är
% koefficienterna för 1, alpha, alpha^2, ... som i poly_fit_ml.
%
% In:
%   fit        anpassningen från poly_fit_ml
%   alphaPrev  senast utvärderade steglängd
%   alphaMax   trust region (0, alphaMax]
%   isBetter   true om alphaPrev gav en förbättring
% Ut:
%   alphaNew   nästa steglängd att utvärdera
%   alphaMax   trust region, krympt om alphaPrev inte gav en förbättring

    if ~isBetter
        alphaMax = alphaPrev/2;                       % krymp trust region
    end
    a = alphaMax*(1:1000)'/1000;                      % rutnät i (0, alphaMax]
    [~, i] = max(a.^(0:numel(fit.y) - 1)*fit.y);      % polynomets värde i rutnätet
    alphaNew = a(i);
end
