function LR = lr_statistic(fit0, fit1, m)
% lr_statistic - Likelihoodkvot för att gå från grad q till q+1 (PROTOTYP).
%
%   LR = m * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaS0^2/sigmaS1^2) )
%
% Tidigare användes faktorn (m - q - 2) som korrektion för få punkter.
%
%   fit0  ML-anpassning av grad q (från poly_fit_ml)
%   fit1  ML-anpassning av grad q+1
%   m     antal punkter
%   LR    stort värde = grad q+1 är klart bättre

    LR = m*(log(fit0.sigmaF^2/fit1.sigmaF^2) + log(fit0.sigmaG^2/fit1.sigmaG^2));
end
