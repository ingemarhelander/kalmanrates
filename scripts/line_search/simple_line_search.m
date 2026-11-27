function [alpha, xNew, fNew, gNew, evals] = simple_line_search(x, f, g, d, K, mu, varInd)
% SIMPLE_LINE_SEARCH - Utvärderar startpunkten och tre punkter per sida längs d.

    alpha = 0;
    xNew = x;
    fNew = f;
    gNew = g;

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.A = 0;
    evals.fit = {};  evals.fitInfo = {};  evals.fitAlpha = {};
    evals.fitF = {};  evals.fitSlope = {};  evals.fitAlphaMax = {};
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';  evals.fits = {};

    gNorm = norm(g);
    dNorm = norm(d);

    alphaSpacing = 1/(gNorm*dNorm);

    alphas = [0, -3:-1, 1:3]*alphaSpacing;
    for k = 2:numel(alphas)
        xTry = x + alphas(k)*d;
        %fprintf('Trying alpha = %.3e, xTry = [%.3e, %.3e, %.3e]\n', alphas(k), xTry);
        [fTry, gTry] = kalmanObj(xTry, K, mu, varInd);
        evals.X(:, end+1) = xTry;
        evals.F(end+1) = fTry;
        evals.G(:, end+1) = gTry;
        evals.Q(end+1) = NaN; % Sammas som i line_search_main.m
        evals.AHat(end+1) = NaN; % Sammas som i line_search_main.m
    end
    evals.A = alphas; % Sammas som i line_search_main.m
end
