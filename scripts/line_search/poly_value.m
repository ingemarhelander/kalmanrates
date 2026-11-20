function p = poly_value(fit, alpha)
% poly_value - Polynomets värde p(alpha).
%
% fit.y är koefficienterna för 1, alpha, alpha^2, ... (från poly_fit_ml).
% Bara poly_value och poly_slope känner till det formatet.
% Utdata har samma form som alpha.

    p = reshape(alpha(:).^(0:numel(fit.y) - 1)*fit.y, size(alpha));
end
