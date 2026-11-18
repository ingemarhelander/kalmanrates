function dp = poly_slope(fit, alpha)
% poly_slope - Derivatan p'(alpha) av det anpassade polynomet längs linjen.
%
% Antar samma format på fit.y som poly_value: koefficienterna för
% 1, alpha, alpha^2, ... som i poly_fit_ml. Ändras formatet behöver båda
% funktionerna ändras.
%
% In:
%   fit    anpassningen från poly_fit_ml
%   alpha  steglängder (godtycklig form)
% Ut:
%   dp     p'(alpha), samma form som alpha

    q = numel(fit.y) - 1;
    if q < 1
        dp = zeros(size(alpha));
        return
    end
    dy = (1:q)'.*fit.y(2:end);          % koefficienter för p', 1, alpha, ...
    dp = reshape(alpha(:).^(0:q - 1)*dy, size(alpha));
end
