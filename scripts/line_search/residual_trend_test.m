function [rho, pValue] = residual_trend_test(alpha, r)
% residual_trend_test - Tests linear trend in residuals.

    alpha = alpha(:);
    r = r(:);

    if numel(alpha) ~= numel(r)
        error('alpha and residuals must have the same length.');
    end

    n = numel(alpha);

    if n < 3
        rho = NaN;
        pValue = NaN;
        return;
    end

    % Correlation between alpha and residuals
    C = corrcoef(alpha, r);
    rho = C(1,2);

    if ~isfinite(rho)
        pValue = NaN;
        return;
    end

    if abs(rho) >= 1
        pValue = 0;
        return;
    end

    % Two-sided Pearson correlation test
    t = rho * sqrt((n-2)/(1-rho^2));
    pValue = betainc((n-2)/(n-2+t^2), (n-2)/2, 0.5);

end