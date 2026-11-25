function alphaNew = find_step_length_newton(fit, nGrid, alphaMax)
% find_step_length_newton - Polynomets första maximum med Newton-Raphson.
%
% Alternativ till find_step_length (ls_maxMethod newton). Ingen kalmanObj
% körs här, bara polynomräkning.
%
%   1. p' räknas i nGrid jämnt fördelade punkter i (0, alphaMax].
%   2. Det första intervallet [lo, hi] där p' går från + till - innehåller
%      det första maximumet.
%   3. Newton-Raphson på p' = 0, t <- t - p'(t)/p''(t), från mittpunkten.
%      Intervallet krymps efter tecknet på p' i varje iterat. Hamnar
%      Newtonsteget utanför intervallet, eller är p'' >= 0, tas
%      mittpunkten i stället (bisektion), så att det alltid konvergerar.
%
% Ändpunkter som i find_step_length:
%   - p' > 0 i alla punkter: alphaMax.
%   - p'(0) <= 0 (polynomet faller redan från start): mittpunkten av
%     första intervallet, alphaMax/(2*nGrid).
%
% In:   fit       ML-anpassningen (från select_poly_degree)
%       nGrid     antal punkter där p' räknas för att hitta intervallet
%       alphaMax  längsta tillåtna steg (1 från line_search_main)
% Ut:   alphaNew  nästa steglängd att köra kalmanObj i

    maxIter = 20;
    tol = 1e-12*alphaMax;

    steps = (1:nGrid)/nGrid*alphaMax;
    k = find(poly_slope(fit, steps) < 0, 1);
    if isempty(k)
        alphaNew = alphaMax;
        return
    end
    if k == 1
        lo = 0;
    else
        lo = steps(k - 1);
    end
    hi = steps(k);
    if poly_slope(fit, lo) <= 0
        % Bara möjligt för lo = 0
        alphaNew = (lo + hi)/2;
        return
    end

    t = (lo + hi)/2;
    for it = 1:maxIter
        s = poly_slope(fit, t);
        if s > 0
            lo = t;
        else
            hi = t;
        end
        c = poly_curvature(fit, t);
        tNew = t - s/c;
        % Skydd: bisektion. Marginalen tol behövs när nollstället ligger
        % precis på en gridpunkt: Newton landar då ett avrundningsfel
        % utanför intervallet och skulle annars bisektera i alla iterationer.
        if ~(c < 0) || tNew < lo - tol || tNew > hi + tol
            tNew = (lo + hi)/2;
        end
        done = abs(tNew - t) <= tol;
        t = tNew;
        if done
            break
        end
    end
    alphaNew = t;
end
