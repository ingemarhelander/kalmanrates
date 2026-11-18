function [q, fit, info] = select_poly_degree(alpha, F, S)
% select_poly_degree - Väljer polynomgrad med likelihoodkvottest (LR).
%
% Börjar med q = 2 och höjer graden så länge grad q+1 är signifikant
% bättre (LR > crit). Högsta grad är m - 2 (minst 2, högst 6), så att
% polynomet aldrig kan gå exakt genom alla F och S.
%
% In:   alpha, F, S  steglängder, funktionsvärden och lutningar g'*d
% Ut:   q     vald grad
%       fit   ML-anpassningen för grad q (från poly_fit_ml)
%       info  LR (ett värde per prövat steg q -> q+1), crit, qMax

    qMin = 2;        % grad 1 har inget maximum
    qCap = 6;
    crit = 3.8415;   % chi2inv(0.95, 1): 5 % risk att höja graden av rent brus

    % Sortera punkterna efter alpha
    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);
    qMax = max(qMin, min(qCap, m - 2));

    q = qMin;
    fit = poly_fit_ml(alpha, F, S, q);
    info = struct('LR', [], 'crit', crit, 'qMax', qMax);
    while q < qMax
        % Grad q+1 går inte att bestämma (rangbrist): behåll q
        try
            fitNext = poly_fit_ml(alpha, F, S, q + 1);
        catch err
            if strcmp(err.identifier, 'poly_fit_ml:rankDeficient')
                break
            end
            rethrow(err);
        end
        LR = lr_statistic(fit, fitNext, m);
        info.LR(end+1) = LR;
        if LR <= crit
            break                    % q+1 är inte signifikant bättre
        end
        q = q + 1;
        fit = fitNext;
    end
end

function LR = lr_statistic(fit0, fit1, m)
% Likelihoodkvot för att gå från grad q (fit0) till q+1 (fit1):
%   LR = (m - q - 2) * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaG0^2/sigmaG1^2) )
% Faktorn (m - q - 2) i stället för m kompenserar för att ML underskattar
% sigma när punkterna är få. Annars höjs graden ofta av brus.
    q = numel(fit0.y) - 1;
    LR = (m - q - 2)*(log(fit0.sigmaF^2/fit1.sigmaF^2) + log(fit0.sigmaG^2/fit1.sigmaG^2));
end
