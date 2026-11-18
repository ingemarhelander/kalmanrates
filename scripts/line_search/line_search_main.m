function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, opts)
% line_search_main - Linjesökning med ML-polynom längs d.
%
% Samma in- och utdata som line_search.m (plus opts), så att
% kalmanOptimizeMain kan välja mellan dem.
%
% Gång:
%   1. Första punkten är alpha = 1, längsta tillåtna steget
%      (find_search_direction skalar d så att alla varianser är positiva).
%   2. En punkt som inte är bättre än F0 och ligger bortom den bästa
%      punkten blir ny övre gräns alphaMax.
%   3. Ingen punkt bättre än F0 än: nästa steg är alphaMax/10.
%   4. Annars: anpassa ML-polynom (select_poly_degree) till punkterna i
%      [0, alphaMax], i t = alpha/alphaMax. find_step_length ger polynomets
%      maximum, som blir nästa punkt.
%   5. Stopp när maximum ligger nära en mätt punkt (< dMin*alphaMax) där
%      |S| <= c*|S0|, när maximum är exakt en mätt punkt, eller efter
%      maxEvals punkter.
% Returnerar den giltiga punkten med störst uppmätt F (alpha = 0 om ingen
% punkt var bättre än startpunkten).
%
% Skriver en tabellrad per punkt: uppmätt F-F0 och dF/dalpha till vänster,
% polynomets grad, nästa alpha och förutsagd ökning till höger.
%
% evals: en kolumn per punkt, startpunkten först.
%   X, F, G  punkt, funktionsvärde, gradient
%   A        steglängd (0 för startpunkten)
%   Q        vald grad (NaN om ingen anpassning gjordes)
%   AHat     nästa steglängd efter punkten
%   reason   varför sökningen stoppade
%
% opts (valfri, ls_-raderna i kalmanOptimizeMain.config):
%   h         (1e-3)  steg längs polynomet i find_step_length, relativt alphaMax
%   dMin      (1e-2)  avstånd relativt alphaMax som räknas som samma punkt
%   c         (0.05)  lutningen räknas som liten när |S| <= c*|S0|
%   maxEvals  (20)    högsta antal kalmanObj-körningar

    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';
    alphas = 0;                  % steglängd för varje kolumn i evals
    alphaMax = 1;                % övre gräns
    alphaNew = alphaMax;         % första punkten
    S0 = g'*d;                   % startlutning
    stop = false;

    fprintf('  Linjesökning: uppmätt startlutning dF/dalpha(0) = %+.4e, första steg alpha = %.4e\n', ...
            S0, alphaNew);
    fprintf('        %-36s | %s\n', 'UPPMÄTT (kalmanObj)', 'POLYNOM (anpassat till punkterna hittills)');
    fprintf('    %2s  %-10s  %-11s  %-11s | %-8s  %-11s  %s\n', ...
            'nr', 'alpha', 'F-F0', 'dF/dalpha', 'polygrad', 'nästa alpha', 'p(nästa)-F0');

    for k = 1:opts.maxEvals
        % Kör Kalman i den nya punkten och spara den
        [fTry, gTry] = kalmanObj(x + alphaNew*d, K, mu, varInd);
        evals.X(:, end+1) = x + alphaNew*d;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        % Radens början skrivs direkt, så att punkten syns även vid krasch
        fprintf('    %2d  %.4e  %+.4e  %+.4e |', k, alphaNew, fTry - f, gTry'*d);

        % Bästa giltiga punkten hittills
        ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
        [fBest, iBest] = max(evals.F(ok));
        okIdx = find(ok);
        alphaBest = alphas(okIdx(iBest));        % 0 om ingen förbättring
        improved = fBest > f;

        % Sämre (eller NaN/Inf) punkt bortom den bästa: ny övre gräns
        if ~(ok(end) && fTry > f) && alphaNew > alphaBest
            alphaMax = min(alphaMax, alphaNew);
        end

        % Ingen förbättring än: krymp steget med faktor 10
        if ~improved
            alphaNew = alphaMax/10;
            evals.AHat(end) = alphaNew;
            if ok(end)
                fprintf(' ingen förbättring än, krymper till alpha = %.4e\n', alphaNew);
            else
                fprintf(' NaN/Inf, krymper till alpha = %.4e\n', alphaNew);
            end
            continue
        end

        % ML-polynom i t = alpha/alphaMax till punkterna i [0, alphaMax].
        % Skalningen behövs när stegen ligger på mycket olika skalor.
        inFit = ok & alphas <= alphaMax*(1 + 1e-12);
        t = alphas(inFit)/alphaMax;
        St = alphaMax*(evals.G(:, inFit)'*d);    % lutning dF/dt
        [q, fit] = select_poly_degree(t, evals.F(inFit), St);
        evals.Q(end) = q;
        tHat = find_step_length(fit, opts.h, 1);
        alphaHat = tHat*alphaMax;
        evals.AHat(end) = alphaHat;
        fprintf(' %8d  %.4e   %+.4e\n', q, alphaHat, poly_value(fit, tHat) - f);

        % Stopp: maximum nära en mätt punkt där lutningen är liten
        [dist, i] = min(abs(alphaHat - alphas));
        slope = evals.G(:, i)'*d;
        if dist < opts.dMin*alphaMax && ok(i) && abs(slope) <= opts.c*abs(S0)
            evals.reason = sprintf(['polynomets maximum alpha = %.4g ligger nära alpha = %.4g, ' ...
                                    'där uppmätt |dF/dalpha| = %.3g <= %.2g * startlutningen'], ...
                                   alphaHat, alphas(i), abs(slope), opts.c);
            stop = true;
            break
        end
        % Stopp: maximum är exakt en mätt punkt. kalmanObj är deterministisk,
        % så en ny körning där ger ingen ny information.
        if dist <= 1e-12*alphaHat
            evals.reason = sprintf(['polynomets maximum alpha = %.4g är redan utvärderat, men uppmätt ' ...
                                    '|dF/dalpha| = %.3g > %.2g * startlutningen'], alphaHat, abs(slope), opts.c);
            stop = true;
            break
        end
        alphaNew = alphaHat;
    end
    if ~stop
        evals.reason = sprintf('maxEvals = %d punkter utvärderade', opts.maxEvals);
    end
    evals.A = alphas;

    % Returnera den giltiga punkten med störst uppmätt F
    ok = find(isfinite(evals.F) & all(isfinite(evals.G), 1));
    [~, i] = max(evals.F(ok));
    j = ok(i);
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);

    % Sammanfattning
    n = numel(alphas) - 1;
    if n == 1
        word = 'utvärdering';
    else
        word = 'utvärderingar';
    end
    fprintf('  Stopp efter %d %s: %s\n', n, word, evals.reason);
    if alpha > 0
        fprintf('  Valt: punkt nr %d, alpha = %.4e, uppmätt F-F0 = %+.4e\n', j - 1, alpha, fNew - f);
    else
        fprintf('  Valt: startpunkten (ingen förbättring)\n');
    end
end

function opts = default_opts(opts)
% Fyller i standardvärden. Okända fält ger fel, så att ett felstavat eller
% gammalt namn i config-filen inte ignoreras tyst.
    def = struct('h', 1e-3, 'dMin', 1e-2, 'c', 0.05, 'maxEvals', 20);
    unknown = setdiff(fieldnames(opts), fieldnames(def));
    if ~isempty(unknown)
        error('Okänd parameter för linjesökningen: %s', strjoin(unknown, ', '));
    end
    defNames = fieldnames(def);
    for i = 1:numel(defNames)
        if ~isfield(opts, defNames{i})
            opts.(defNames{i}) = def.(defNames{i});
        end
    end
    if ~(isscalar(opts.h) && isnumeric(opts.h) && opts.h > 0)
        error('ls_h måste vara ett positivt tal.');
    end
end
