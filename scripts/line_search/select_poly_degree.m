function [q, fit, info] = select_poly_degree(alpha, F, S, qPrev)
% select_poly_degree - Väljer polynomgrad längs linjen (Req. 25, 26, 28-30).
%
% Provar q = 2, 3, ... och väljer den lägsta grad där residualerna ser ut som
% brus. En för låg grad lämnar ett slätt mönster, vilket ger stort positivt rho.
% Bara ett för stort rho höjer graden (ensidigt test).
%
%   alpha  steglängder (bara punkter där F och S är ändliga)
%   F      funktionsvärden
%   S      lutning längs linjen, g'*d
%   qPrev  grad från förra varvet, eller [] första gången (valfri)
%   q      vald grad, fit är anpassningen för den graden
%   info   (valfri) en rad per provad grad: info.q, info.rho = [rho_F rho_S],
%          info.tau, info.fits. Används för felsökning och figurer. tau är
%          NaN för sista graden om info inte efterfrågas.

    qMin   = 2;     % lägsta grad, används när punkterna är få
    qCap   = 6;     % högsta grad
    mMin   = 4;     % färre punkter än så: alltid qMin
    level  = 0.95;  % bara brus överskrider tau i 5 % av fallen
    nSim   = 1000;  % antal brusförsök för att ta fram tau
    margin = 0.1;   % en lägre grad än förra varvet måste vara klart bättre (hysteres)

    if nargin < 4 || isempty(qPrev)
        qPrev = 0;  % ingen grad räknas som lägre
    end

    % rho jämför grannar, så punkterna måste ligga i ordning
    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);

    % Minst två fler observationer (2m) än koefficienter (q+1)
    qMax = qMin;
    if m >= mMin
        qMax = min(qCap, 2*m - 3);
    end

    tolExact = 1e-10;   % residualer under tolExact*data räknas som avrundningsfel
    data = [F - mean(F), S];
    info = struct('q', [], 'rho', zeros(0, 2), 'tau', [], 'fits', {{}});

    for q = qMin:qMax
        fit = poly_fit_ml(alpha, F, S, q, false);   % false: inga figurer
        r = [fit.rF(:), fit.rG(:)];
        % Exakt anpassning: residualerna är bara avrundningsfel
        isExact = norm(r(:)) <= tolExact*norm(data(:));
        % En grupp (F eller S) med bara avrundningsfel kvar räknas som brus.
        % Avrundningsfel är inte vitt brus och kan annars ge stort rho.
        rho = lag1_autocorr(r);
        rho(vecnorm(r) <= tolExact*vecnorm(data)) = 0;
        % Stanna vid högsta graden, vid exakt anpassning eller när rho är litet.
        % Under förra varvets grad krävs att rho är minst margin under tau.
        isLast = q == qMax || isExact;
        tau = NaN;
        if ~isLast || nargout > 2
            tau = whiteness_threshold(alpha, q, fit, level, nSim);
        end
        info.q(end+1) = q;  info.rho(end+1, :) = rho;
        info.tau(end+1) = tau;  info.fits{end+1} = fit;
        if isLast || max(rho) <= tau - margin*(q < qPrev)
            break
        end
    end
end

function tau = whiteness_threshold(alpha, q, fit, level, nSim)
% Gräns för rho: anpassa samma polynom till rent brus många gånger och se hur
% stort rho blir. Bara brus överskrider tau i 5 % av fallen.
% Fast frö ger samma tau varje gång och påverkar inte rng.

    % Designmatris för polynom av grad q. alpha skalas till [-1, 1] för
    % noggrannhetens skull, det ändrar inte vilka polynom som kan anpassas.
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
