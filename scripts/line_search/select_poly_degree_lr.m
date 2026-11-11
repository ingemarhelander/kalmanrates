function [q, fit, info] = select_poly_degree_lr(alpha, F, S)
% select_poly_degree_lr - PROTOTYP: väljer polynomgrad med likelihoodkvottest.
%
% Alternativ till select_poly_degree (autokorrelation) för få punkter.
% Provar q = 2, 3, ... och höjer graden så länge q+1 ger en signifikant
% bättre ML-anpassning än q:
%
%   LR = (m - q - 2) * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaS0^2/sigmaS1^2) )
%
% där index 0 är grad q och 1 är grad q+1. Utan faktorn (m - q - 2), alltså
% med m, är det det vanliga LR = 2(l1 - l0). Faktorn är en småstickprovs-
% korrektion: ML underskattar sigma med få frihetsgrader, och utan den höjs
% graden för ofta (ca 30 % falsklarm med 5-6 punkter i simulering).
% Graden höjs om LR > 3.84 (95 %-kvantilen för chi2 med 1 frihetsgrad).
%
% Högsta grad är m - 2, så att varken F eller S kan anpassas exakt (då går
% sigma mot 0 och LR mot oändligheten).
%
%   alpha, F, S  som i select_poly_degree
%   q, fit       vald grad och ML-anpassningen för den
%   info         info.q (provade grader), info.LR (LR för q -> q+1),
%                info.crit, info.fits

    crit = 3.8415;   % chi2inv(0.95, 1), hårdkodad för att slippa Statistics Toolbox
    qMin = 2;
    qCap = 6;

    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);
    qMax = max(qMin, min(qCap, m - 2));

    q = qMin;
    fit = poly_fit_ml(alpha, F, S, q, false);
    info = struct('q', q, 'LR', [], 'crit', crit, 'fits', {{fit}});
    while q < qMax
        fitNext = poly_fit_ml(alpha, F, S, q + 1, false);
        LR = (m - q - 2)*(log(fit.sigmaF^2/fitNext.sigmaF^2) + log(fit.sigmaG^2/fitNext.sigmaG^2));
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
