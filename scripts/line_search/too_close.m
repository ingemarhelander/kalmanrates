function tooClose = too_close(alphaHat, alphas, dMinRel)
% too_close - Ligger alphaHat för nära en redan utvärderad punkt (Req. 31)?
%
% Gemensam definition för evaluate (som stoppar sökningen när en förbättring
% finns) och line_search_main (som annars provar ett kortare steg).
%
% alphaHat är för nära om
%   1. avståndet till ett utvärderat steg är mindre än dMinRel*alphaHat.
%      Relativt det föreslagna steget och inte det största, som kan vara ett
%      misslyckat långt steg: då hindras sökningen från att förfina i den
%      skala där förbättringen finns,
%   2. alphaHat är mindre än dMinRel gånger det minsta utvärderade steget.
%      Då ligger det i praktiken i startpunkten alpha = 0. Avståndet dit är
%      alphaHat självt, så det fallet fångas aldrig av 1.
%
%   alphaHat  föreslagen nästa steglängd
%   alphas    steglängd för alla utvärderade punkter, startpunkten (0) och
%             ogiltiga punkter (NaN/Inf) inräknade
%   dMinRel   opts.dMinRel. Mindre värde ger färre punkter som är för nära.

    steps = alphas(alphas > 0);
    tooClose = min(abs(alphaHat - steps)) < dMinRel*alphaHat ...
               || alphaHat < dMinRel*min(steps);
end
