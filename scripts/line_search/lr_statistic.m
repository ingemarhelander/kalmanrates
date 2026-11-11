function LR = lr_statistic(fit0, fit1, m, q)
% lr_statistic - Likelihoodkvot för att gå från grad q till q+1 (PROTOTYP).
%
%   LR = (m - q - 2) * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaS0^2/sigmaS1^2) )
%
% Faktorn (m - q - 2) i stället för m är en korrektion för få punkter,
% eftersom ML underskattar sigma när frihetsgraderna är få.
%
%   fit0  ML-anpassning av grad q (från poly_fit_ml)
%   fit1  ML-anpassning av grad q+1
%   m     antal punkter
%   q     grad för fit0
%   LR    stort värde = grad q+1 är klart bättre

    LR = (m - q - 2)*(log(fit0.sigmaF^2/fit1.sigmaF^2) + log(fit0.sigmaG^2/fit1.sigmaG^2));
end
