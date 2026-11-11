function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, method, lsOpts)
% line_search_main - Linjesökning med polynom (designspec 5.2.4).
%
% Ska ersätta line_search.m när den är klar, och har därför samma in- och utdata.
% Varje varv: utvärdera en ny punkt, anpassa ML-polynomet (ml_poly_main) och
% hitta nästa steglängd alphaHat. Stoppar när evaluate säger till
% (komponent 10) eller när maxEvals punkter är utvärderade.
% Returnerar den giltiga utvärderade punkt som har störst uppmätt F.
%
% Skriver ut en rad per utvärderad punkt.
%
% method (valfri): gradval i ml_poly_main, 'lag1' (standard) eller 'lr'.
% lsOpts (valfri): struct som ersätter standardvärden i opts till evaluate,
%   t.ex. lsOpts.minPoints = 6 eller lsOpts.c = 0.1.
%
% evals: en kolumn per utvärderad punkt, startpunkten först (X, F, G).
%   evals.A       steglängd för varje punkt (startpunkten har 0)
%   evals.Q       vald grad efter varje punkt (NaN om ingen anpassning gjordes)
%   evals.AHat    föreslagen nästa steglängd efter varje punkt (NaN om ingen)
%   evals.reason  varför sökningen stoppade

    if nargin < 8 || isempty(method)
        method = 'lag1';
    end
    maxEvals = 20;     % högsta antal nya punkter per sökriktning
    opts.dMinRel = 1e-3;  % minsta avstånd till en utvärderad punkt, relativt största steget (Req. 31)
    opts.c       = 0.5;   % stopp när |f'_j| <= c |f'_0| (lutningen har planat ut)
    if nargin >= 9 && ~isempty(lsOpts)
        for name = fieldnames(lsOpts)'
            opts.(name{1}) = lsOpts.(name{1});
        end
    end

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';
    alphas = 0;                              % steglängd för varje kolumn i evals
    alphaNew = initial_step_length(g, d);      % första nya punkten (Req. 22)
    alphaMax = 1;                            % trust region (0, alphaMax]
    qPrev = [];                              % grad från förra varvet (hysteres)
    stop = false;                            % true när evaluate stoppat sökningen

    for k = 1:maxEvals
        xTry = x + alphaNew*d;
        [fTry, gTry] = kalmanObj(xTry, K, mu, varInd);
        evals.X(:, end+1) = xTry;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
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
        [qPrev, fit] = ml_poly_main(alphas(ok), evals.F(ok), evals.G(:, ok)'*d, qPrev, method);
        evals.Q(end) = qPrev;

        % isBetter styr bara trust region i find_step_length, inte när vi stoppar
        isBetter = isfinite(fTry) && fTry > f;
        [alphaHat, alphaMax] = find_step_length(fit, alphaNew, alphaMax, isBetter);
        evals.AHat(end) = alphaHat;

        % Stoppkriterier (komponent 10), prövas innan alphaHat utvärderas
        [stop, reason] = evaluate(fit, alphaHat, alphas, ok, evals.F, gTry'*d, g'*d, opts);
        if stop
            fprintf('    line search stop: %s\n', reason);
            evals.reason = reason;
            break
        end
        alphaNew = alphaHat;
    end
    if ~stop
        evals.reason = sprintf('maxEvals = %d punkter utvärderade', maxEvals);
        fprintf('    line search stop: %s\n', evals.reason);
    end
    evals.A = alphas;

    % Returnera den giltiga punkt som har störst uppmätt F. Startpunkten är
    % med, så utan förbättring returneras den.
    ok = find(isfinite(evals.F) & all(isfinite(evals.G), 1));
    [~, i] = max(evals.F(ok));
    j = ok(i);
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);
end
