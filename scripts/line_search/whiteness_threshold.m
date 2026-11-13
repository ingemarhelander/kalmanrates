function tau = whiteness_threshold(alpha, q, fit, level, nSim)
% whiteness_threshold - Gräns tau för lag-1-autokorrelationen (Req. 29).
%
% Anpassar samma polynom till rent brus många gånger och ser hur stort rho
% blir. Bara brus överskrider tau i (1 - level) av fallen.
% Fast frö ger samma tau varje gång och påverkar inte rng.
%
%   alpha  steglängder, sorterade
%   q      polynomgrad
%   fit    anpassningen från poly_fit_ml (vikterna fit.wF och fit.wG används)
%   level  kvantil, t.ex. 0.95 (valfri)
%   nSim   antal brusförsök, t.ex. 1000 (valfri)

    if nargin < 4, level = 0.95; end
    if nargin < 5, nSim = 1000; end

    % Designmatris för polynom av grad q. alpha skalas till [-1, 1] för
    % noggrannhetens skull, det ändrar inte vilka polynom som kan anpassas.
    alpha = alpha(:);
    m = numel(alpha);
    s = (alpha(end) - alpha(1))/2;
    x = (alpha - alpha(1))/s - 1;
    k = 1:q;
    A = [x.^(0:q); zeros(m, 1), k.*x.^(k - 1)/s];               % rader för F, sedan S

    Aw = sqrt([fit.wF*ones(m, 1); fit.wG*ones(m, 1)]).*A;       % vikta raderna så att F och S väger lika
    E  = randn(RandStream('mt19937ar', 'Seed', 0), 2*m, nSim);  % rent brus, en kolumn per försök
    R  = E - Aw*(Aw\E);                                         % residualer när datan är rent brus
    rho = sort(max(lag1_autocorr(R(1:m, :)), lag1_autocorr(R(m+1:end, :))));
    tau = rho(ceil(level*nSim));
end
