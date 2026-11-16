function [q, fit, info] = select_poly_degree_lr(alpha, F, S)
% select_poly_degree_lr - Väljer polynomgrad med likelihoodkvottest (PROTOTYP).
%
% Alternativ till select_poly_degree_lag1. Provar q = 2, 3, ... och höjer
% graden så länge grad q+1 är signifikant bättre (lr_statistic > crit).
% Högsta grad är m - 2, så att varken F eller S kan anpassas exakt.
% Ingen hysteres (qPrev).
%
%   alpha  steglängder (bara punkter där F och S är ändliga)
%   F      funktionsvärden
%   S      lutning längs linjen, g'*d
%   q      vald grad, fit är ML-anpassningen för den graden
%   info   (valfri) info.q (godkända grader i ordning), info.LR (LR för
%          varje prövat steg q -> q+1), info.crit, info.qMax, info.fits

    qMin = 2;        % lägsta grad
    qCap = 6;        % högsta grad
    crit = 3.8415;   % chi2inv(0.95, 1): bara brus ger LR > crit i 5 % av fallen

    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);
    qMax = max(qMin, min(qCap, m - 2));

    q = qMin;
    fit = poly_fit_ml(alpha, F, S, q, false);
    info = struct('q', q, 'LR', [], 'crit', crit, 'qMax', qMax, 'fits', {{fit}});
    while q < qMax
        fitNext = poly_fit_ml(alpha, F, S, q + 1, false);
        LR = lr_statistic(fit, fitNext, m);
        info.LR(end+1) = LR;
        if LR <= crit
            break                    % q+1 är inte signifikant bättre
        end
        q = q + 1;
        fit = fitNext;
        info.q(end+1) = q;
        info.fits{end+1} = fit;
    end
end
