function [q, fit, info] = select_poly_degree(alpha, F, S, opts)
% select_poly_degree - Väljer polynomgrad med likelihoodkvottest (LR).
%
% Börjar med q = 2 och höjer graden så länge grad q+1 är signifikant
% bättre (LR > crit). Är q+1 inte signifikant bättre prövas även q+2 mot q
% med en strängare gräns (crit2). En term kan nämligen sakna effekt medan
% nästa behövs, t.ex. den udda grad 3-termen för en funktion som är nästan
% symmetrisk kring intervallets mitt. Högsta grad är m - 2 (minst 2, högst
% 6), så att polynomet aldrig kan gå exakt genom alla F och S.
%
% Hoppet q -> q+2 prövas bara när m >= q + 5. Med färre punkter är
% korrektionsfaktorn (m - q - 2) i LR för liten och LR följer chi2 dåligt:
% med m = 6 höjdes graden av rent brus i 13 % av fallen i stället för 5 %.
%
% In:   alpha, F, S  steglängder, funktionsvärden och lutningar g'*d
%       opts         (valfri) struct med
%                    lookAhead (true)  false: pröva bara q -> q+1, som
%                                      före hoppet q -> q+2 infördes
% Ut:   q     vald grad
%       fit   ML-anpassningen för grad q (från poly_fit_ml)
%       info  crit, crit2, qMax, m, och ett värde per prövat steg
%             qFrom -> qTo: LR, qFrom, qTo, critTest (gränsen som
%             användes), sigmaF och sigmaG (rader: [qFrom; qTo])

    qMin = 2;        % grad 1 har inget maximum
    qCap = 6;
    crit = 3.8415;   % chi2inv(0.95, 1): 5 % risk att höja graden av rent brus
    crit2 = 9.2103;  % chi2inv(0.99, 2): hopp två grader, strängare eftersom
                     % ett extra test annars ger fler falsklarm
    lookAhead = true;
    if nargin >= 4 && isfield(opts, 'lookAhead')
        lookAhead = logical(opts.lookAhead);
    end

    % Sortera punkterna efter alpha
    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);
    qMax = max(qMin, min(qCap, m - 2));

    q = qMin;
    fit = poly_fit_ml(alpha, F, S, q);
    info = struct('LR', [], 'qFrom', [], 'qTo', [], 'critTest', [], ...
                  'sigmaF', zeros(2, 0), 'sigmaG', zeros(2, 0), ...
                  'crit', crit, 'crit2', crit2, 'qMax', qMax, 'm', m);
    while q < qMax
        % Varmstart i vikterna för grad q, så att LR >= 0
        [fitNext, ok] = try_fit(alpha, F, S, q + 1, fit);
        if ~ok
            break                    % grad q+1 går inte att bestämma: behåll q
        end
        [info, raise] = lr_test(info, fit, fitNext, m, crit);
        if raise
            q = q + 1;
            fit = fitNext;
            continue
        end

        % q+1 är inte signifikant bättre: pröva q+2 mot q, om det finns
        % tillräckligt många punkter
        if ~lookAhead || q + 2 > qMax || m < q + 5
            break
        end
        [fitNext2, ok] = try_fit(alpha, F, S, q + 2, fitNext);
        if ~ok
            break
        end
        [info, raise] = lr_test(info, fit, fitNext2, m, crit2);
        if ~raise
            break
        end
        q = q + 2;
        fit = fitNext2;
    end
end

function [fit, ok] = try_fit(alpha, F, S, q, fit0)
% poly_fit_ml med varmstart i fit0. ok = false vid rangbrist.
    fit = [];
    ok = true;
    try
        fit = poly_fit_ml(alpha, F, S, q, fit0);
    catch err
        if ~strcmp(err.identifier, 'poly_fit_ml:rankDeficient')
            rethrow(err);
        end
        ok = false;
    end
end

function [info, raise] = lr_test(info, fit0, fit1, m, critTest)
% LR-test från fit0 till fit1 mot gränsen critTest. Sparar testet i info.
    LR = lr_statistic(fit0, fit1, m);
    info.LR(end+1) = LR;
    info.qFrom(end+1) = numel(fit0.y) - 1;
    info.qTo(end+1) = numel(fit1.y) - 1;
    info.critTest(end+1) = critTest;
    info.sigmaF(:, end+1) = [fit0.sigmaF; fit1.sigmaF];
    info.sigmaG(:, end+1) = [fit0.sigmaG; fit1.sigmaG];
    raise = LR > critTest;
end

function LR = lr_statistic(fit0, fit1, m)
% Likelihoodkvot för att gå från grad q (fit0) till en högre grad (fit1):
%   LR = (m - q - 2) * ( log(sigmaF0^2/sigmaF1^2) + log(sigmaG0^2/sigmaG1^2) )
% Faktorn (m - q - 2) i stället för m kompenserar för att ML underskattar
% sigma när punkterna är få. Annars höjs graden ofta av brus.
    q = numel(fit0.y) - 1;
    LR = (m - q - 2)*(log(fit0.sigmaF^2/fit1.sigmaF^2) + log(fit0.sigmaG^2/fit1.sigmaG^2));
end
