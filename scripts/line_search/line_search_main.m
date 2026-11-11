function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd)
% line_search_main - Linjesökning med polynom (designspec 5.2.4).
%
% Ska ersätta line_search.m när den är klar, och har därför samma in- och utdata.
% Varje varv: utvärdera en ny punkt, anpassa ML-polynomet (ml_poly_main) och
% hitta nästa steglängd alphaHat. Stoppar om evaluate säger till eller när
% maxEvals punkter är utvärderade.
% Returnerar den utvärderade punkt som har störst anpassat värde p(alpha_j).
%
% Skriver ut en rad per utvärderad punkt.
%
% evals: en kolumn per utvärderad punkt, startpunkten först (X, F, G),
% och steglängden för varje punkt i evals.A (startpunkten har 0).

    maxEvals = 20;     % högsta antal nya punkter per sökriktning

    evals.X = x;  evals.F = f;  evals.G = g;
    alphas = 0;                              % steglängd för varje kolumn i evals
    alphaNew = initial_step_length(g, d);      % första nya punkten (Req. 22)
    alphaMax = 1;                            % trust region (0, alphaMax]
    qPrev = [];                              % grad från förra varvet (hysteres)
    fit = [];                                % anpassning till alla giltiga punkter

    for k = 1:maxEvals
        xTry = x + alphaNew*d;
        [fTry, gTry] = kalmanObj(xTry, K, mu, varInd);
        evals.X(:, end+1) = xTry;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        fprintf('    line search eval %2d: alpha = %.4e, f-f0 = %+.4e, f'' = %+.4e\n', ...
                k, alphaNew, fTry - f, gTry'*d);

        % Punkter med NaN eller Inf tas inte med i anpassningen (Req. 42).
        % En sådan punkt är aldrig en förbättring: halvera steget och krymp
        % trust region, så att polynomsteget inte går tillbaka dit.
        ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
        if ~ok(end) || sum(ok) < 2
            alphaNew = alphaNew/2;
            alphaMax = min(alphaMax, alphaNew);
            continue
        end
        [qPrev, fit] = ml_poly_main(alphas(ok), evals.F(ok), evals.G(:, ok)'*d, qPrev);

        % isBetter styr bara trust region i find_step_length, inte när vi stoppar
        isBetter = isfinite(fTry) && fTry > f;
        [alphaHat, alphaMax] = find_step_length(fit, alphaNew, alphaMax, isBetter);

        % Stoppkriterier, prövas innan alphaHat utvärderas
        if evaluate(fit, alphaHat, alphas, ok)
            break
        end
        alphaNew = alphaHat;
    end

    % Returnera den giltiga punkt som har störst anpassat värde p(alpha_j),
    % som påverkas mindre av brus än F. fit bygger alltid på alla giltiga
    % punkter. Utan anpassning (färre än två giltiga punkter): startpunkten.
    j = 1;
    if ~isempty(fit)
        ok = find(isfinite(evals.F) & all(isfinite(evals.G), 1));
        [~, i] = max(poly_value(fit, alphas(ok)));
        j = ok(i);
    end
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);
end
