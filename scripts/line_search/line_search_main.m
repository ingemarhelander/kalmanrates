function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, opts)
% line_search_main - Linjesökning med ML-polynom.
%
% Har samma in- och utdata som line_search.m (plus det valfria argumentet
% opts), så att kalmanOptimizeMain kan välja mellan dem.
%
% Sökningen stoppar aldrig utan att ha hittat en punkt med högre F än
% startpunkten (utom när maxEvals tar slut). d är en uppåtriktning
% (dF/dalpha(0) > 0), så ett tillräckligt kort steg ger alltid en förbättring.
%
%   1. Första punkten är det längsta tillåtna steget, alpha = 1.
%      find_search_direction skalar d så att alpha = 1 håller alla
%      varianser positiva.
%   2. En punkt som inte är bättre än startpunkten (eller är NaN/Inf) och
%      ligger bortom den bästa punkten hittills blir ny övre gräns alphaMax.
%   3. Så länge ingen punkt är bättre än startpunkten: nästa steg är
%      alphaMax/10. Med en extremt brant start (t.ex. en varians nära noll)
%      kan maximum ligga många tiopotenser under 1, och då behövs det.
%   4. När en bättre punkt finns tar polynomet över. ML-polynomet anpassas
%      till de giltiga punkterna i [0, alphaMax], i den skalade variabeln
%      t = alpha/alphaMax. Punkterna bortom alphaMax är alla sämre än
%      startpunkten, och utan skalningen blir anpassningen numeriskt omöjlig
%      när stegen ligger på mycket olika skalor. find_step_length går längs
%      polynomet från t = 0 med steg h tills p' blir negativ första gången
%      och föreslår mittpunkten mellan de två sista stegen.
%   5. Stopp om förslaget ligger närmare än dMin*alphaMax från en utvärderad
%      punkt OCH den uppmätta lutningen där är liten, |S| <= c |S_0|
%      (S = g'*d). Stopp också om förslaget är exakt en redan utvärderad
%      punkt: kalmanObj är deterministisk, så en ny körning där ger ingen
%      ny information. Annars körs kalmanObj i förslaget och vi går
%      tillbaka till 2.
%
% Returnerar den giltiga utvärderade punkt som har störst uppmätt F.
% Startpunkten är med, så om maxEvals tar slut utan förbättring returneras
% alpha = 0.
%
% Skriver ut en tabell med en rad per utvärderad punkt. Vänstra delen är
% uppmätt med kalmanObj (F-F0 och riktningsderivatan dF/dalpha = g'*d),
% högra delen kommer från polynomet (grad, föreslagen nästa steglängd och
% förutsagd ökning där), eller anger att steget krymps.
%
% evals: en kolumn per utvärderad punkt, startpunkten först (X, F, G).
%   evals.A       steglängd för varje punkt (startpunkten har 0)
%   evals.Q       vald grad efter varje punkt (NaN om ingen anpassning gjordes)
%   evals.AHat    nästa steglängd efter varje punkt (NaN för startpunkten)
%   evals.reason  varför sökningen stoppade
%
% opts (valfri, sätts i kalmanOptimizeMain.config med prefixet ls_). Fält
% som saknas får standardvärdet:
%   h            (1e-3)  steglängd när find_step_length går längs polynomet,
%                        relativt alphaMax
%   dMin         (1e-2)  avstånd till en utvärderad punkt, relativt alphaMax,
%                        som räknas som samma punkt i stoppvillkoret
%   c            (0.05)  lutningen räknas som liten när |S| <= c |S_0|
%   maxEvals     (20)    högsta antal nya punkter per sökriktning
%   degreeMethod ('lr')  gradval i ml_poly_main: 'lr', 'score',
%                        'threshold' eller 'lag1'

    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';
    alphas = 0;                  % steglängd för varje kolumn i evals
    alphaMax = 1;                % övre gräns: längsta tillåtna steg
    alphaNew = alphaMax;         % första punkten: längsta tillåtna steg
    qPrev = [];                  % grad från förra varvet (hysteres i 'lag1')
    S0 = g'*d;                   % startlutning
    stop = false;

    fprintf('  Linjesökning: uppmätt startlutning dF/dalpha(0) = %+.4e, första steg alpha = %.4e\n', ...
            S0, alphaNew);
    fprintf('        %-36s | %s\n', 'UPPMÄTT (kalmanObj)', 'POLYNOM (anpassat till punkterna hittills)');
    fprintf('    %2s  %-10s  %-11s  %-11s | %-8s  %-11s  %s\n', ...
            'nr', 'alpha', 'F-F0', 'dF/dalpha', 'polygrad', 'nästa alpha', 'p(nästa)-F0');

    for k = 1:opts.maxEvals
        [fTry, gTry] = kalmanObj(x + alphaNew*d, K, mu, varInd);
        evals.X(:, end+1) = x + alphaNew*d;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        % Raden avslutas när nästa steg är känt. Början skrivs redan här så
        % att punkten syns även om anpassningen kraschar.
        fprintf('    %2d  %.4e  %+.4e  %+.4e |', k, alphaNew, fTry - f, gTry'*d);

        ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
        [fBest, iBest] = max(evals.F(ok));
        okIdx = find(ok);
        alphaBest = alphas(okIdx(iBest));        % 0 om ingen förbättring
        improved = fBest > f;

        % Punkt utan förbättring bortom den bästa punkten: ny övre gräns
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

        % Polynomet i t = alpha/alphaMax, anpassat till punkterna i [0, alphaMax]
        inFit = ok & alphas <= alphaMax*(1 + 1e-12);
        t = alphas(inFit)/alphaMax;
        St = alphaMax*(evals.G(:, inFit)'*d);    % lutning dF/dt
        [qPrev, fit] = ml_poly_main(t, evals.F(inFit), St, qPrev, opts.degreeMethod);
        evals.Q(end) = qPrev;
        tHat = find_step_length(fit, opts.h, 1);
        alphaHat = tHat*alphaMax;
        evals.AHat(end) = alphaHat;
        fprintf(' %8d  %.4e   %+.4e\n', qPrev, alphaHat, poly_value(fit, tHat) - f);

        % Stopp: polynomets maximum ligger vid en punkt vi redan har mätt,
        % och den uppmätta lutningen där är liten
        [dist, i] = min(abs(alphaHat - alphas));
        slope = evals.G(:, i)'*d;
        if dist < opts.dMin*alphaMax && ok(i) && abs(slope) <= opts.c*abs(S0)
            evals.reason = sprintf(['polynomets maximum alpha = %.4g ligger nära alpha = %.4g, ' ...
                                    'där uppmätt |dF/dalpha| = %.3g <= %.2g * startlutningen'], ...
                                   alphaHat, alphas(i), abs(slope), opts.c);
            stop = true;
            break
        end
        % Förslaget är exakt en punkt vi redan har mätt (men lutningen där är
        % inte liten). kalmanObj är deterministisk, så en ny körning där ger
        % samma värden och polynomet ändras inte.
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

    % Returnera den giltiga punkt som har störst uppmätt F. Startpunkten är
    % med, så utan förbättring returneras den.
    ok = find(isfinite(evals.F) & all(isfinite(evals.G), 1));
    [~, i] = max(evals.F(ok));
    j = ok(i);
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);

    % Sammanfattning: varför sökningen stoppade och vilken punkt som valdes
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
% Fyller i standardvärden för fält som saknas. Okända fält ger fel, så att
% ett felstavat (eller gammalt) namn i config-filen inte ignoreras i tysthet.
    def = struct('h', 1e-3, 'dMin', 1e-2, 'c', 0.05, 'maxEvals', 20, 'degreeMethod', 'lr');
    names = fieldnames(opts);
    unknown = setdiff(names, fieldnames(def));   % fungerar även när opts saknar fält
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
