function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, opts)
% line_search_main - Linjesökning med ML-polynom längs d.
%
% Samma in- och utdata som line_search.m (plus opts), så att
% kalmanOptimizeMain kan välja mellan dem.
%
% Gång:
%   0. dF/dalpha(0) <= 0: d är ingen uppåtriktning, startpunkten returneras
%      utan någon kalmanObj-körning.
%   1. alphaMax = min(opts.alphaMax, 0.99 * steget till första undre gränsen
%      längs d), om opts.xLbd är given. Längre steg gör en varians negativ,
%      och komplext F, för alla större alpha (x + alpha*d är linjärt).
%      Första punkten ges av initial_step_length, i [alphaMin, alphaMax].
%      Komplext F (t.ex. utan xLbd) sparas som NaN.
%   2. En punkt som inte är bättre än F0 (eller är ogiltig) och ligger
%      bortom den bästa punkten blir ny övre gräns alphaMax.
%   3. Ingen punkt bättre än F0 än: nästa steg är alphaMax/10.
%   4. Annars: anpassa ML-polynom (select_poly_degree) till de giltiga
%      punkterna i [0, alphaMax], i t = alpha/alphaMax. find_step_length ger
%      polynomets maximum, som blir nästa punkt.
%   5. Är maximum exakt en mätt punkt som är ogiltig eller sämre än den
%      bästa: gå halvvägs mellan den bästa punkten och närmaste mätta
%      granne på den sida som lutningen i den bästa punkten pekar mot.
%      Ligger de redan närmare än dMin*alphaMax: stopp.
%   6. Stopp när lutningen i den bästa punkten är liten, |S| <= c*|S0|,
%      eller polynomets maximum ligger nära en mätt punkt med liten lutning
%      (båda tidigast efter minEvals punkter), när maximum är den bästa
%      punkten, eller efter maxEvals punkter.
% Returnerar den giltiga punkten med störst uppmätt F (alpha = 0 om ingen
% punkt var bättre än startpunkten).
%
% Skriver en tabellrad per punkt: alpha, uppmätt F-F0 och dF/dalpha (eller
% komplex/NaN/Inf), polynomets grad, nästa alpha och varför ('krymp /10',
% 'polynomets max', 'halvvägs', 'stopp'). Sist varför sökningen stoppade
% och vilken punkt som valdes.
%
% evals: en kolumn per giltig punkt, startpunkten först. Ogiltiga punkter
% (NaN/Inf/komplex) tas bort innan evals returneras, eftersom denoise
% använder evals.
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
%   minEvals  (1)     stoppet på liten lutning gäller först efter så många punkter
%   alphaMax  (20)    längsta tillåtna steg
%   alphaMin  (0.1)   kortaste första steg
%   xLbd      ([])    undre gränser för x (sätts av kalmanOptimizeMain, inte
%                     i config). Tom: ingen gränskontroll

    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';
    alphas = 0;                  % steglängd för varje kolumn i evals
    alphaMax = opts.alphaMax;    % övre gräns: längsta tillåtna steg
    alphaBound = Inf;            % 0.99 * steget till första undre gränsen
    if ~isempty(opts.xLbd)
        alphaBound = max_feasible_step(x, d, opts.xLbd);
        alphaMax = min(alphaMax, alphaBound);
    end
    alphaNew = initial_step_length(g, d, alphaMax, opts.alphaMin);   % första punkten
    S0 = g'*d;                   % startlutning
    stop = false;

    if alphaBound < opts.alphaMax
        boundNote = ' (variansgräns)';
    else
        boundNote = '';
    end
    fprintf('  Linjesökning: dF/dalpha(0) = %+.3e, alphaMax = %.3g%s\n', S0, alphaMax, boundNote);
    fprintf('    %2s  %-9s  %-10s  %-10s  %4s  %-9s  %s\n', ...
            'nr', 'alpha', 'F-F0', 'dF/dalpha', 'grad', 'nästa', 'varför');

    % d är ingen uppåtriktning: inget steg längs d kan förväntas ge högre F,
    % så ingen Kalman-körning görs och startpunkten returneras
    nLoop = opts.maxEvals;
    if ~(S0 > 0)
        nLoop = 0;
        stop = true;
        evals.reason = 'dF/dalpha(0) <= 0, d är ingen uppåtriktning';
    elseif ~(alphaMax > 0)
        nLoop = 0;
        stop = true;
        evals.reason = 'en varians ligger redan på sin gräns och d pekar nedåt (alphaMax = 0)';
    end

    for k = 1:nLoop
        [fTry, gTry] = kalmanObj(x + alphaNew*d, K, mu, varInd);
        % Komplext F (en varians har blivit negativ, log i barriärtermen)
        % sparas som NaN, så att evals.F förblir reell
        isComplex = any(imag(fTry) ~= 0) || any(imag(gTry) ~= 0);
        if isComplex
            fTry = NaN;  gTry = NaN(size(gTry));
        else
            fTry = real(fTry);  gTry = real(gTry);
        end
        evals.X(:, end+1) = x + alphaNew*d;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        % Raden avslutas när nästa steg är känt (print_rest). Början skrivs
        % redan här så att punkten syns även om anpassningen kraschar.
        if isComplex
            fprintf('    %2d  %.3e  %-10s  %-10s', k, alphaNew, 'komplex', '');
        elseif ~(isfinite(fTry) && all(isfinite(gTry)))
            fprintf('    %2d  %.3e  %-10s  %-10s', k, alphaNew, 'NaN/Inf', '');
        else
            fprintf('    %2d  %.3e  %+.3e  %+.3e', k, alphaNew, fTry - f, gTry'*d);
        end

        ok = valid_points(evals);
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
            print_rest([], alphaNew, 'krymp /10');
            continue
        end

        % Stopp: den uppmätta lutningen i den bästa punkten är liten
        % (först när minst minEvals nya punkter är utvärderade)
        sBest = evals.G(:, okIdx(iBest))'*d;
        if k >= opts.minEvals && abs(sBest) <= opts.c*abs(S0)
            evals.reason = sprintf('lutningen i bästa punkten alpha = %.3g är liten, |dF/dalpha| <= %.2g*|dF/dalpha(0)|', ...
                                   alphaBest, opts.c);
            print_rest([], [], 'stopp');
            stop = true;
            break
        end

        % ML-polynom i t = alpha/alphaMax till de giltiga punkterna i
        % [0, alphaMax] (även de med sämre F än F0). Punkter längre ut tas
        % inte med: de ligger ofta tiopotenser bort, och ett polynom kan inte
        % följa både dem och toppen nära 0. Skalningen behövs när stegen
        % ligger på mycket olika skalor.
        inFit = ok & alphas <= alphaMax*(1 + 1e-12);
        t = alphas(inFit)/alphaMax;
        St = alphaMax*(evals.G(:, inFit)'*d);    % lutning dF/dt
        [q, fit] = select_poly_degree(t, evals.F(inFit), St);
        evals.Q(end) = q;
        tHat = find_step_length(fit, opts.h, 1);
        alphaHat = tHat*alphaMax;
        evals.AHat(end) = alphaHat;

        % Stopp: polynomets maximum ligger vid en punkt vi redan har mätt,
        % och den uppmätta lutningen där är liten
        [dist, i] = min(abs(alphaHat - alphas));
        slope = evals.G(:, i)'*d;
        % (först när minst minEvals nya punkter är utvärderade)
        if k >= opts.minEvals && dist < opts.dMin*alphaMax && ok(i) && abs(slope) <= opts.c*abs(S0)
            evals.reason = sprintf('polynomets max ligger nära alpha = %.3g, där |dF/dalpha| <= %.2g*|dF/dalpha(0)|', ...
                                   alphas(i), opts.c);
            print_rest(q, [], 'stopp');
            stop = true;
            break
        end
        % Förslaget är exakt en punkt vi redan har mätt, och den punkten är
        % ogiltig (NaN/komplex) eller sämre än den bästa. Polynomet ser inte
        % den punkten (den ligger på gränsen alphaMax eller är ogiltig), så
        % gå halvvägs mellan den bästa punkten och närmaste utvärderade
        % punkt på den sida där den uppmätta lutningen i den bästa punkten
        % säger att maximum ligger (positiv: höger, annars vänster). Den nya
        % punkten ligger då mellan två mätta grannar och är aldrig redan
        % utvärderad.
        if dist <= 1e-12*alphaHat && i ~= okIdx(iBest)
            if sBest > 0
                neighbour = min(alphas(alphas > alphaBest));
                if isempty(neighbour)
                    neighbour = alphaMax;    % ingen mätt punkt till höger
                end
            else
                neighbour = max(alphas(alphas < alphaBest));
            end
            % Stopp: den bästa punkten och grannen ligger redan närmare än
            % dMin*alphaMax (t.ex. när F ökar ända fram till en variansgräns)
            if abs(neighbour - alphaBest) < opts.dMin*alphaMax
                evals.reason = sprintf(['maximum inringat mellan alpha = %.4g och %.4g ' ...
                                        '(avstånd < dMin*alphaMax)'], alphaBest, neighbour);
                print_rest(q, [], 'stopp');
                stop = true;
                break
            end
            alphaNew = (alphaBest + neighbour)/2;
            evals.AHat(end) = alphaNew;
            print_rest(q, alphaNew, 'halvvägs');
            continue
        end
        % Förslaget är exakt den bästa punkten (men lutningen där är inte
        % liten). kalmanObj är deterministisk, så en ny körning där ger
        % samma värden och polynomet ändras inte.
        if dist <= 1e-12*alphaHat
            evals.reason = sprintf('polynomets max är den bästa punkten alpha = %.3g, redan utvärderad', alphaHat);
            print_rest(q, [], 'stopp');
            stop = true;
            break
        end
        alphaNew = alphaHat;
        print_rest(q, alphaNew, 'polynomets max');
    end
    if ~stop
        evals.reason = sprintf('maxEvals = %d punkter utvärderade', opts.maxEvals);
    end
    evals.A = alphas;

    % Returnera den giltiga punkt som har störst uppmätt F. Startpunkten är
    % med, så utan förbättring returneras den.
    ok = find(valid_points(evals));
    [~, i] = max(evals.F(ok));
    j = ok(i);
    alpha = alphas(j);
    xNew = evals.X(:, j);
    fNew = evals.F(j);
    gNew = evals.G(:, j);
    nInvalid = numel(alphas) - numel(ok);

    % Ta bort ogiltiga punkter (NaN/Inf/komplex) ur evals innan den lämnas
    % ut. Loopen behövde dem för att krympa alphaMax, men denoise använder
    % evals och klarar inte NaN.
    evals.X = evals.X(:, ok);  evals.F = evals.F(ok);  evals.G = evals.G(:, ok);
    evals.A = evals.A(ok);  evals.Q = evals.Q(ok);  evals.AHat = evals.AHat(ok);

    % Sammanfattning: varför sökningen stoppade och vilken punkt som valdes
    n = numel(alphas) - 1;
    if n == 1
        word = 'utvärdering';
    else
        word = 'utvärderingar';
    end
    invalidNote = '';
    if nInvalid > 0
        invalidNote = sprintf(', %d ogiltiga', nInvalid);
    end
    if alpha > 0
        chosen = sprintf('valt nr %d, alpha = %.3e, F-F0 = %+.3e', j - 1, alpha, fNew - f);
    else
        chosen = 'valt startpunkten (ingen förbättring)';
    end
    fprintf('  Stopp efter %d %s%s: %s\n  -> %s\n', n, word, invalidNote, evals.reason, chosen);
end

function print_rest(q, alphaNext, why)
% Avslutar tabellraden: polynomets grad, nästa alpha och varför.
% Tomma q eller alphaNext skrivs som blanka.
    qStr = '';
    if ~isempty(q)
        qStr = sprintf('%d', q);
    end
    nextStr = '';
    if ~isempty(alphaNext)
        nextStr = sprintf('%.3e', alphaNext);
    end
    fprintf('  %4s  %-9s  %s\n', qStr, nextStr, why);
end

function alphaLbd = max_feasible_step(x, d, xLbd)
% 0.99 gånger längsta steget längs d innan någon variabel når sin undre
% gräns (samma marginal som i find_search_direction). Inf om ingen gräns nås.
    hits = (d < 0) & isfinite(xLbd);
    alphaLbd = 0.99*min([Inf; -(x(hits) - xLbd(hits))./d(hits)]);
end

function ok = valid_points(evals)
% Punkter som får användas: ändliga F och G (komplexa är redan NaN).
    ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
end

function opts = default_opts(opts)
% Fyller i standardvärden. Okända fält ger fel, så att ett felstavat eller
% gammalt namn i config-filen inte ignoreras tyst.
    def = struct('h', 1e-3, 'dMin', 1e-2, 'c', 0.05, 'maxEvals', 20, ...
                 'alphaMax', 20, 'alphaMin', 0.1, 'minEvals', 1, 'xLbd', []);
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
