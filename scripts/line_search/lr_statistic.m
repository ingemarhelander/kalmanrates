function LR = lr_statistic(fit0, fit1, m)
% lr_statistic - Likelihoodkvot för att gå från grad q till q+1 (PROTOTYP).
%
%   LR = (m - q - 2) * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaS0^2/sigmaS1^2) )
%
% Faktorn (m - q - 2) i stället för m är en korrektion för få punkter,
% eftersom ML underskattar sigma när frihetsgraderna är få. Med bara m
% höjs graden ofta av brus när punkterna är få (t.ex. m = 6, q = 3 -> 4
% ger faktor 6 i stället för 1), se testSelectPolyDegreeLr/cubicGivesThree.
%
%   fit0  ML-anpassning av grad q (från poly_fit_ml), q = numel(fit0.y) - 1
%   fit1  ML-anpassning av grad q+1
%   m     antal punkter
%   LR    stort värde = grad q+1 är klart bättre

    q = numel(fit0.y) - 1;
    LR = (m - q - 2)*(log(fit0.sigmaF^2/fit1.sigmaF^2) + log(fit0.sigmaG^2/fit1.sigmaG^2));
end
