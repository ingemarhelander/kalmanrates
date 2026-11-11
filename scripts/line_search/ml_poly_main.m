function [q, fit, info] = ml_poly_main(alpha, F, S, qPrev)
% ml_poly_main - ML-polynom längs linjen (ML-boxen i line search).
%
% Ingång till ML-delen. Anropar select_poly_degree, som anpassar ett
% ML-polynom (poly_fit_ml) för varje grad och väljer graden med
% autokorrelationen (lag1_autocorr).
%
%   alpha  steglängder (bara punkter där F och S är ändliga)
%   F      funktionsvärden
%   S      lutning längs linjen, g'*d
%   qPrev  grad från förra varvet, eller [] första gången (valfri)
%   q      vald grad, fit är ML-anpassningen för den graden
%   info   (valfri) detaljer per provad grad, se select_poly_degree

    if nargin < 4
        qPrev = [];
    end
    if nargout > 2
        [q, fit, info] = select_poly_degree(alpha, F, S, qPrev);
    else
        [q, fit] = select_poly_degree(alpha, F, S, qPrev);
    end
end
