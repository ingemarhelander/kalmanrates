function rho = lag1_autocorr(r)
% lag1_autocorr - Lag-1-autokorrelation för residualer (Req. 27).
%
% Mäter om residualer som ligger bredvid varandra brukar ha samma tecken.
%   rho nära 0  bara brus är kvar
%   rho nära 1  ett slätt mönster är kvar, alltså är graden för låg
%
% r ska vara sorterad efter alpha. Varje kolumn är en egen följd och ger ett
% eget värde i rho.

    if isrow(r)
        r = r(:);                  % konvertera till kolumnvektor
    end
    r = r - mean(r, 1);            % centrera hela kolumnerna, så att rho inte påverkas av konstant bias
    rho = sum(r(1:end-1, :).*r(2:end, :), 1) ./ sum(r.^2, 1);
    rho(~isfinite(rho)) = 0;       % alla residualer noll: räknas som brus
end
