function fit = poly_fit_ml(alpha, F, S, q, fit0)
% poly_fit_ml - ML-anpassning av ett polynom av grad q till både
% funktionsvärden F och lutningar S (ekv. 34, 38, 39 i ML-dokumentet).
%
% Bruset på F och S skattas separat. Iterationen:
%   1. Starta med vikterna wF = wG = 1, eller med vikterna i fit0.
%   2. Skatta koefficienterna med viktad minsta kvadrat på F och S.
%   3. Skatta sigmaF^2 och sigmaG^2 från residualerna, sätt w = 1/sigma^2.
%   4. Upprepa 2-3 tills vikterna slutar ändras.
%
% In:   alpha, F, S  steglängder, funktionsvärden och lutningar g'*d (m x 1)
%       q            polynomgrad
%       fit0         (valfri) anpassning av lägre grad till samma punkter.
%                    Iterationen startar i dess vikter. Varje varv höjer
%                    likelihooden, så resultatet blir minst lika bra som
%                    fit0 (LR >= 0) och hamnar i samma läge som fit0 i
%                    stället för att följa F eller S beroende på graden.
% Ut:   fit.y              koefficienter för 1, alpha, alpha^2, ... (q+1 x 1)
%       fit.rF, fit.rG     residualer för F och S
%       fit.sigmaF/sigmaG  skattade standardavvikelser
%       fit.wF, fit.wG     vikter 1/sigma^2
%       fit.iterations, fit.converged
%
% Ger felet 'poly_fit_ml:rankDeficient' om punkterna inte bestämmer
% alla koefficienter.

    alpha = alpha(:);
    F = F(:);
    S = S(:);
    m = length(alpha);

    if length(F) ~= m || length(S) ~= m
        error('All input vectors must have the same number of observations.');
    end
    if q < 1 || q ~= floor(q)
        error('q must be a positive integer.');
    end

    % Anpassa i t = alpha/alphaScale, annars blir alpha^q för litet vid
    % korta steg. Lutningen i t är S*alphaScale.
    alphaScale = max(abs(alpha));
    if alphaScale == 0
        alphaScale = 1;
    end
    S = S*alphaScale;

    % Designmatriser: a*y = p(t), b*y = p'(t)
    a = bsxfun(@power, alpha/alphaScale, 0:q);
    b = [zeros(m,1), ...
         bsxfun(@times, bsxfun(@power, alpha/alphaScale, 0:q-1), 1:q)];

    if rank([a; b]) < q+1
        error('poly_fit_ml:rankDeficient', ['The design matrix is rank deficient. ', ...
               'The selected step lengths cannot identify all coefficients.']);
    end

    wF = 1;
    wG = 1;
    if nargin >= 5
        % fit0.wG gäller lutningen i alpha; här är lutningen i t
        wF = fit0.wF;
        wG = fit0.wG/alphaScale^2;
    end
    maxIter = 1000;
    tol = 1e-8;
    converged = false;

    % Golv för varianserna (avrundningsfel), annars kan vikterna hoppa
    % mellan varven när residualerna är nästan noll.
    floorF = (1e-10*max(sqrt(mean(F.^2)), 1))^2;
    floorG = (1e-10*max(sqrt(mean(S.^2)), 1))^2;

    for k = 1:maxIter
        % Viktad minsta kvadrat på F och S tillsammans
        y = [sqrt(wF)*a; sqrt(wG)*b] \ [sqrt(wF)*F; sqrt(wG)*S];

        % Nya varianser och vikter från residualerna
        rF = F - a*y;
        rG = S - b*y;
        sigmaF2 = max(mean(rF.^2), floorF);
        sigmaG2 = max(mean(rG.^2), floorG);
        wF_new = 1/sigmaF2;
        wG_new = 1/sigmaG2;

        change = max(abs([wF_new/wF - 1, wG_new/wG - 1]));
        wF = wF_new;
        wG = wG_new;
        if change < tol
            converged = true;
            break
        end
    end

    if ~converged
        fprintf('    poly_fit_ml: grad %d konvergerade inte på %d iterationer\n', q, maxIter);
    end

    % Sista y i loopen skattades med föregående vikter: skatta om med de
    % slutliga
    y = [sqrt(wF)*a; sqrt(wG)*b] \ [sqrt(wF)*F; sqrt(wG)*S];
    rF = F - a*y;
    rG = S - b*y;
    sigmaF2 = max(mean(rF.^2), floorF);
    sigmaG2 = max(mean(rG.^2), floorG);

    % Räkna om från t till alpha
    fit.y = y ./ alphaScale.^(0:q)';
    fit.rF = rF;
    fit.rG = rG/alphaScale;
    fit.sigmaF = sqrt(sigmaF2);
    fit.sigmaG = sqrt(sigmaG2)/alphaScale;
    fit.wF = 1/sigmaF2;
    fit.wG = 1/fit.sigmaG^2;
    fit.iterations = k;
    fit.converged = converged;

    fit.alphaScale = alphaScale;
    fit.yScaled = y;
end
