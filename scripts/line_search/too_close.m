function tf = too_close(alphaHat, alphas, dMinRel)
% too_close - Ligger alphaHat för nära en redan utvärderad punkt (Req. 31)?
%
% Avståndet mäts relativt största utvärderade steget, max(alphas).
% alphas innehåller alla utvärderade steglängder (startpunkten 0 inräknad).

    tf = min(abs(alphaHat - alphas)) < dMinRel*max(alphas);
end
