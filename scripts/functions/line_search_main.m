function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd)
% line_search_main - Linjesökning med polynom (designspec 5.2.4).
%
% Ska ersätta line_search.m när den är klar, och har därför samma in- och utdata.
% Varje varv: utvärdera en ny punkt, välj polynomgrad, anpassa polynomet och
% hitta nästa steglängd. Klart när ett polynomsteg ger en förbättring.
%
% evals: en kolumn per utvärderad punkt, startpunkten först (X, F, G).

    maxEvals = 20;   % högsta antal nya punkter

    evals.X = x;  evals.F = f;  evals.G = g;
    alphas = 0;                              % steglängd för varje kolumn i evals
    alphaNew = first_step_length(g, d);      % första nya punkten (Req. 22)
    alphaMax = 1;                            % trust region (0, alphaMax]
    qPrev = [];                              % grad från förra varvet (hysteres)

    for k = 1:maxEvals
        xTry = x + alphaNew*d;
        [fTry, gTry] = kalmanObj(xTry, K, mu, varInd);
        evals.X(:, end+1) = xTry;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        alphas(end+1) = alphaNew; %#ok<AGROW>

        % Första punkten godkänns inte direkt, annars blir det backtracking
        isBetter = isfinite(fTry) && fTry > f;
        if k > 1 && isBetter
            break
        end

        % Punkter med NaN eller Inf tas inte med i anpassningen
        ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
        if sum(ok) < 2
            alphaNew = alphaNew/2;           % bara startpunkten är giltig: korta steget
            continue
        end
        [qPrev, fit] = select_poly_degree(alphas(ok), evals.F(ok), evals.G(:, ok)'*d, qPrev);
        [alphaNew, alphaMax] = find_step_length(fit, alphaNew, alphaMax, isBetter);
    end

    % Returnera den bästa giltiga punkten, startpunkten om ingen är bättre
    F = evals.F;
    F(~(isfinite(F) & all(isfinite(evals.G), 1))) = -Inf;
    [~, j] = max(F);
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);
end
