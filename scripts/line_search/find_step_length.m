function [alphaNew, alphaMax] = find_step_length(fit, alphaPrev, alphaMax, isBetter)
% find_step_length - Nästa steglängd från ML-polynomet.
%
% Söker efter den första lokala maximipunkten inom trust region
% [0, alphaMax]. Om ingen sådan punkt hittas väljs alphaMax.
%
% In:
%   fit        Anpassningen från poly_fit_ml
%   alphaPrev  Senast utvärderade steglängd
%   alphaMax   Aktuell övre gräns för trust region
%   isBetter   true om alphaPrev gav en förbättring
%
% Ut:
%   alphaNew   Nästa steglängd att utvärdera
%   alphaMax   Uppdaterad trust region

    %% Uppdatera trust region

    if ~isBetter
        alphaMax = alphaPrev / 2;
    end

    %% Hämta polynomkoefficienter

    y = fit.y(:);
    q = numel(y) - 1;

    if q < 1
        error('The fitted polynomial must have degree at least 1.');
    end

    %% Beräkna första- och andraderivatan

    % P'(alpha)
    dy = (1:q)' .* y(2:end);

    % P''(alpha)
    if q >= 2
        ddy = (2:q)' .* (1:q-1)' .* y(3:end);
    else
        ddy = [];
    end

    %% Hitta stationära punkter

    % Nollställen till P'(alpha)
    stationaryPoints = roots(flipud(dy));

    % Behåll endast reella rötter
    tol = 1e-10 * max(1, alphaMax);

    isReal = abs(imag(stationaryPoints)) <= tol;
    stationaryPoints = real(stationaryPoints(isReal));

    % Behåll endast punkter inom trust region
    insideRegion = stationaryPoints >= -tol & ...
                   stationaryPoints <= alphaMax + tol;

    stationaryPoints = stationaryPoints(insideRegion);

    % Korrigera små numeriska avvikelser vid intervallets gränser
    stationaryPoints = min(max(stationaryPoints, 0), alphaMax);

    % Sortera i stigande ordning
    stationaryPoints = sort(stationaryPoints);

    %% Hitta första lokala maximum

    alphaNew = alphaMax;

    if q >= 2 && ~isempty(stationaryPoints)

        % Beräkna andraderivatan i alla stationära punkter
        secondDerivative = polyval(flipud(ddy), stationaryPoints);

        % Ett lokalt maximum kräver P''(alpha) < 0
        isMaximum = secondDerivative < 0;

        % Välj den första stationära punkten som är ett lokalt maximum
        firstMaximum = find(isMaximum, 1, 'first');

        if ~isempty(firstMaximum)
            alphaNew = stationaryPoints(firstMaximum);
        end
    end

end
