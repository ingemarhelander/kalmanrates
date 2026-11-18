function dp = poly_slope(fit, alpha)
% poly_slope - Polynomets derivata p'(alpha).
%
% Samma format på fit.y som i poly_value. Utdata har samma form som alpha.

    q = numel(fit.y) - 1;
    if q < 1
        dp = zeros(size(alpha));     % konstant polynom
        return
    end
    dy = (1:q)'.*fit.y(2:end);       % koefficienterna för p'
    dp = reshape(alpha(:).^(0:q - 1)*dy, size(alpha));
end
