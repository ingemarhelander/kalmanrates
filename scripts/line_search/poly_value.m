function p = poly_value(fit, alpha)
% poly_value - Värdet av det anpassade polynomet p(alpha) längs linjen.
%
% Enda stället som känner till formatet på fit.y. När poly_fit_ml byts ut
% mot ML-versionen behöver bara den här funktionen ändras.
%
% TILLFÄLLIG VERSION: antar att fit.y är koefficienterna för
% 1, alpha, alpha^2, ... som i den tillfälliga poly_fit_ml.
%
% In:
%   fit    anpassningen från poly_fit_ml
%   alpha  steglängder (godtycklig form)
% Ut:
%   p      p(alpha), samma form som alpha

    p = reshape(alpha(:).^(0:numel(fit.y) - 1)*fit.y, size(alpha));
end
