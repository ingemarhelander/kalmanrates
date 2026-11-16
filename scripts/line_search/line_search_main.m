function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, opts)
% line_search_main - Linjesökning med polynom (designspec 5.2.4).
%
% Ska ersätta line_search.m när den är klar, och har därför samma in- och utdata
% (plus det valfria argumentet opts).
% Varje varv: utvärdera en ny punkt, anpassa ML-polynomet (ml_poly_main) och
% hitta nästa steglängd alphaHat. Stoppar när evaluate säger till
% (komponent 10) eller när maxEvals punkter är utvärderade.
% Returnerar den giltiga utvärderade punkt som har störst uppmätt F.
%
% Skriver ut den uppmätta startlutningen, en tabell med en rad per utvärderad
% punkt och sist varför sökningen stoppade och vilken punkt som valdes.
% Tabellens vänstra del är uppmätt med kalmanObj (F-F0 och riktningsderivatan
% dF/dalpha = g'*d), den högra kommer från det anpassade polynomet (grad,
% föreslagen nästa steglängd och förutsagd ökning där).
%
% evals: en kolumn per utvärderad punkt, startpunkten först (X, F, G).
%   evals.A       steglängd för varje punkt (startpunkten har 0)
%   evals.Q       vald grad efter varje punkt (NaN om ingen anpassning gjordes)
%   evals.AHat    föreslagen nästa steglängd efter varje punkt (NaN om ingen)
%   evals.reason  varför sökningen stoppade
%
% opts (valfri, sätts i kalmanOptimizeMain.config med prefixet ls_). Fält
% som saknas får standardvärdet:
%   maxEvals     (20)   högsta antal nya punkter per sökriktning
%   c            (0.2)  stopp när |f'_j| <= c |f'_0| (lutningen har planat ut)
%   dMinRel      (1e-3) minsta avstånd till en utvärderad punkt, relativt
%                       det föreslagna steget alphaHat (Req. 31)
%   noiseFactor  (0)    stopp när förutsagd ökning < noiseFactor*sigmaF.
%                       0 stänger av kriteriet (tillfälligt standardvärde)
%   minPoints    (0)    stopp på lutning eller brus först när minst så många
%                       punkter (startpunkten inräknad) är utvärderade
%   degreeMethod ('lr') gradval i ml_poly_main, 'lr' eller 'lag1'

    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);
    maxEvals = opts.maxEvals;
    method = opts.degreeMethod;

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';
    alphas = 0;                              % steglängd för varje kolumn i evals
    % trust region (0, alphaMax]. Högst 1: find_search_direction skalar d så
    % att alpha = 1 är det längsta steg som håller alla varianser positiva.
    alphaMax = 1;
    alphaNew = initial_step_length(g, d, alphaMax);      % första nya punkten (Req. 22)
    qPrev = [];                              % grad från förra varvet (hysteres)
    stop = false;                            % true när evaluate stoppat sökningen

    % Startlutningen f'_0 är det lutningsvillkoret jämför mot. Negativ
    % f'_0 betyder att d inte är en ökningsriktning.
    fprintf('  Linjesökning: uppmätt startlutning dF/dalpha(0) = %+.4e, första steg alpha = %.4e\n', ...
            g'*d, alphaNew);
    % Kolumnerna till vänster om | är uppmätta med kalmanObj i punkten alpha.
    % Kolumnerna till höger kommer från polynomet p som anpassats till alla
    % punkter hittills: dess grad, dess förslag på nästa alpha och den
    % ökning det förutsäger där.
    fprintf('        %-36s | %s\n', 'UPPMÄTT (kalmanObj)', 'POLYNOM (anpassat till punkterna hittills)');
    fprintf('    %2s  %-10s  %-11s  %-11s | %-8s  %-11s  %s\n', ...
            'nr', 'alpha', 'F-F0', 'dF/dalpha', 'polygrad', 'nästa alpha', 'p(nästa)-F0');

    for k = 1:maxEvals
        xTry = x + alphaNew*d;
        [fTry, gTry] = kalmanObj(xTry, K, mu, varInd);
        evals.X(:, end+1) = xTry;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
        % Raden avslutas längre ner, när grad och nästa steg är kända. Början
        % skrivs redan här så att punkten syns även om anpassningen kraschar.
        fprintf('    %2d  %.4e  %+.4e  %+.4e |', k, alphaNew, fTry - f, gTry'*d);

        % Punkter med NaN eller Inf tas inte med i anpassningen (Req. 42).
        % En sådan punkt är aldrig en förbättring: halvera steget och krymp
        % trust region, så att polynomsteget inte går tillbaka dit.
        ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
        if ~ok(end) || sum(ok) < 2
            alphaNew = alphaNew/2;
            alphaMax = min(alphaMax, alphaNew);
            fprintf(' ingen anpassning (NaN/Inf), halverar till alpha = %.4e\n', alphaNew);
            continue
        end
        [qPrev, fit] = ml_poly_main(alphas(ok), evals.F(ok), evals.G(:, ok)'*d, qPrev, method);
        evals.Q(end) = qPrev;

        % isBetter styr bara trust region i find_step_length, inte när vi stoppar
        isBetter = isfinite(fTry) && fTry > f;
        [alphaHat, alphaMax] = find_step_length(fit, alphaNew, alphaMax, isBetter);
        % Ingen förbättring än och alphaHat för nära en utvärderad punkt
        % (samma definition som i evaluate, se too_close): prova ett 10 gånger
        % kortare steg i stället för att stoppa med alpha = 0
        note = '';
        if ~(max(evals.F(ok)) > f) && too_close(alphaHat, alphas, opts.dMinRel)
            alphaHat = min(alphas(alphas > 0))/10;
            note = '  (polynomets förslag för nära en utvärderad punkt: minsta alpha/10)';
            if alphaHat < 1e-8*max(alphas)   % ge upp, steget är försumbart
                fprintf(' %8d\n', qPrev);
                evals.reason = 'ingen förbättring, steget har krympts under 1e-8 gånger största steget';
                stop = true;
                break
            end
        end
        evals.AHat(end) = alphaHat;
        fprintf(' %8d  %.4e   %+.4e%s\n', qPrev, alphaHat, poly_value(fit, alphaHat) - f, note);

        % Stoppkriterier (komponent 10), prövas innan alphaHat utvärderas
        [stop, reason] = evaluate(fit, alphaHat, alphas, ok, evals.F, gTry'*d, g'*d, opts);
        if stop
            % Polynomets max ligger på trust region-gränsen, som redan är
            % utvärderad: sökningen har nått taket, och det är därför nästa
            % punkt är "för nära". Säg det i klartext.
            [dist, i] = min(abs(alphaHat - alphas));
            if contains(reason, 'för nära') && alphaHat >= alphaMax*(1 - 1e-9) ...
                    && dist <= opts.dMinRel*alphaHat && ok(i)
                if evals.G(:, i)'*d > 0
                    slope = 'positiv';
                else
                    slope = 'negativ';
                end
                reason = sprintf('nådde taket alpha = %.3g med %s lutning (Req. 31)', alphaMax, slope);
            end
            evals.reason = reason;
            break
        end
        alphaNew = alphaHat;
    end
    if ~stop
        evals.reason = sprintf('maxEvals = %d punkter utvärderade', maxEvals);
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
% ett felstavat namn i config-filen inte ignoreras i tysthet.
    def = struct('maxEvals', 20, 'c', 0.2, 'dMinRel', 1e-3, 'noiseFactor', 0, ...
                 'minPoints', 0, 'degreeMethod', 'lr');
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
end
