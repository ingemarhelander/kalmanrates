function d2p = poly_curvature(fit, alpha)
% poly_curvature - Polynomets andraderivata p''(alpha).
%
% Samma format på fit.y som i poly_value. Utdata har samma form som alpha.

    q = numel(fit.y) - 1;
    if q < 2
        d2p = zeros(size(alpha));    % högst linjärt polynom
        return
    end
    ddy = ((2:q).*(1:q - 1))'.*fit.y(3:end);     % koefficienterna för p''
    d2p = reshape(alpha(:).^(0:q - 2)*ddy, size(alpha));
end
