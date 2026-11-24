function [alpha, xNew, fNew, gNew, evals] = line_search_main(x, f, g, d, K, mu, varInd, opts, isSimple)
% line_search_main - Linjesökning med ML-polynom längs d.
%
% Samma in- och utdata som line_search.m (plus opts), så att
% kalmanOptimizeMain kan välja mellan dem.
%
% Gång:
%   0. Är dF/dalpha(0) <= 0 pekar d inte uppåt: returnera startpunkten.
%   1. alpha ligger alltid i (0, 1]. Riktningen d är redan skalad så att
%      alpha = 1 inte går förbi någon variansgräns (find_search_direction
%      och denoise), så här finns ingen gränskontroll. Första punkten
%      alpha1 från initial_step_length_main, i [alphaMin, 1].
%   1b. Förvalda punkter (initial_points_main): alpha1, 3*alpha1/4,
%      alpha1/2 och alpha1/4 mäts innan något polynom anpassas, så att
%      första anpassningen har 5 punkter (med startpunkten). Då är
%      ML-skattningen av bruset inte degenererad och gradvalet kan välja
%      grad 3 direkt (se initial_points_main).
%      Punkterna skrivs ut som vanliga tabellrader med "förvald punkt".
%   2. Ogiltig punkt (NaN/Inf), eller ingen förbättring, bortom den
%      bästa punkten: ny övre gräns alphaMax (börjar på 1; om
%      shrinkAlphaMax = 1).
%   3. Därefter: ML-polynom (select_poly_degree) genom alla giltiga punkter,
%      även de bortom alphaMax; nästa punkt är polynomets maximum i
%      [0, alphaMax] (find_step_length).
%   4. Stopp när
%        - lutningen i bästa punkten är liten, |S| <= c*|S0| (tidigast
%          efter minEvals punkter; standard 1, alltså direkt),
%        - maximum ligger närmare än dMin*alphaMax en redan mätt punkt,
%        - maxEvals punkter är utvärderade.
% Returnerar den giltiga punkten med störst F (alpha = 0 om ingen var
% bättre än startpunkten).
%
% Utskrift (bara om verbose = 1): en tabellrad per punkt (alpha, F-F0, dF/dalpha, grad, nästa
% alpha, varför), sist varför sökningen stoppade och vilken punkt som valdes.
%
% evals: en kolumn per giltig punkt, startpunkten först. Ogiltiga punkter
% (NaN/Inf) tas bort innan evals returneras, eftersom denoise
% använder evals.
%   X, F, G  punkt, funktionsvärde, gradient
%   A        steglängd (0 för startpunkten)
%   Q        vald grad (NaN om ingen anpassning gjordes)
%   AHat     nästa steglängd efter punkten
%   reason   varför sökningen stoppade
%   fits     (cell, tom om inte opts.saveFits) indata och resultat för
%            varje anrop av select_poly_degree, för replayDegreeSelection:
%            k (punktnummer), t = alpha/alphaMax, F, St (dF/dt), alphaMax,
%            q (vald grad), tHat (polynomets max i t)
%
% opts (valfri, ls_-raderna i kalmanOptimizeMain.config):
%   h         (1e-3)  steg längs polynomet i find_step_length, relativt alphaMax
%   dMin      (1e-2)  avstånd relativt alphaMax som räknas som samma punkt
%   c         (0.05)  lutningen räknas som liten när |S| <= c*|S0|
%   maxEvals  (20)    högsta antal kalmanObj-körningar
%   minEvals  (1)     stoppet på liten lutning gäller först efter så många
%                     punkter. 1 påverkar ingenting (gäller direkt).
%   alphaMin  (0.1)   kortaste första steg
%   showDegree (0)    1: skriv gradvalet under varje tabellrad (LR för varje
%                     prövat steg q -> q+1 eller q -> q+2, sigmaF och sigmaG
%                     före och efter)
%   saveFits  (0)     1: spara indata till varje gradval i evals.fits
%   shrinkAlphaMax (1) 1: en punkt bortom bästa som inte är bättre än F0
%                     blir ny övre gräns alphaMax (maximum söks bara under
%                     den; t, St och dMin skalas med den). 0: alphaMax = 1
%                     hela sökningen.
%   verbose   (0)     0: ingen utskrift (kalmanOptimizeMain skriver en rad per
%                     iteration). 1: tabellen och sammanfattningen nedan.

    

    if nargin < 8 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    if isSimple
        [alpha, xNew, fNew, gNew, evals] = simple_line_search(x, f, g, d, K, mu, varInd);
        return
    end

    evals.X = x;  evals.F = f;  evals.G = g;
    evals.fit = {};  evals.fitInfo = {};  evals.fitAlpha = {};
    evals.fitF = {};  evals.fitSlope = {};  evals.fitAlphaMax = {};
    evals.Q = NaN;  evals.AHat = NaN;  evals.reason = '';  evals.fits = {};
    alphas = 0;                  % steglängd för varje kolumn i evals
    alphaMax = 1;                % övre gräns, krymper när en punkt bortom bästa är dålig
    alpha1 = initial_step_length_main(g, d, opts.alphaMin);     % första punkten, i [alphaMin, 1]
    pending = initial_points_main(alpha1);   % förvalda punkter: alpha1, 3/4, 1/2, 1/4 av alpha1
    S0 = g'*d;                   % startlutning
    stop = false;

    v = opts.verbose;            % 0: ingen utskrift alls (main skriver)
    if v
        fprintf('  Linjesökning: dF/dalpha(0) = %+.3e\n', S0);
        fprintf('    %2s  %-9s  %-10s  %-10s  %4s  %-9s  %s\n', ...
                'nr', 'alpha', 'F-F0', 'dF/dalpha', 'grad', 'nästa', 'varför');
    end

    % Inget steg längs d kan ge högre F: hoppa över loopen och returnera
    % startpunkten
    nLoop = opts.maxEvals;
    if ~(S0 > 0)
        nLoop = 0;
        stop = true;
        evals.reason = 'dF/dalpha(0) <= 0, d är ingen uppåtriktning';
    end

    for k = 1:nLoop
        % De förvalda punkterna mäts först, sedan bestämmer polynomet
        if ~isempty(pending)
            alphaNew = pending(1);
            pending(1) = [];
        end
        [fTry, gTry] = kalmanObj(x + alphaNew*d, K, mu, varInd);
        evals.X(:, end+1) = x + alphaNew*d;  evals.F(end+1) = fTry;  evals.G(:, end+1) = gTry;
        evals.Q(end+1) = NaN;  evals.AHat(end+1) = NaN;
        alphas(end+1) = alphaNew; %#ok<AGROW>
        % Första halvan av tabellraden skrivs nu, så att punkten syns även om
        % anpassningen kraschar; print_rest avslutar raden.
        if ~v
            % ingen tabell
        elseif ~(isfinite(fTry) && all(isfinite(gTry)))
            fprintf('    %2d  %.3e  %-10s  %-10s', k, alphaNew, 'NaN/Inf', '');
        else
            fprintf('    %2d  %.3e  %+.3e  %+.3e', k, alphaNew, fTry - f, gTry'*d);
        end

        % Bästa giltiga punkten hittills (startpunkten om ingen är bättre)
        ok = valid_points(evals);
        Fvalid = evals.F;
        Fvalid(~ok) = -Inf;                      % ogiltiga kan aldrig bli bäst
        [~, iBest] = max(Fvalid);                % iBest: kolumn i evals
        alphaBest = alphas(iBest);

        % Ogiltig punkt, eller ingen förbättring, bortom bästa punkten:
        % ny övre gräns (om shrinkAlphaMax = 1; med 0 ligger alphaMax kvar på 1)
        if opts.shrinkAlphaMax && ~(ok(end) && fTry > f) && alphaNew > alphaBest
            alphaMax = min(alphaMax, alphaNew);
        end

        % Förvalda punkter kvar: mät nästa innan något polynom anpassas
        if ~isempty(pending)
            evals.AHat(end) = pending(1);
            print_rest(v, [], pending(1), 'förvald punkt');
            continue
        end

        % ML-polynom i t = alpha/alphaMax genom alla giltiga punkter (utan
        % NaN/Inf), även de sämre än F0 och de bortom alphaMax (t > 1).
        % Maximum söks ändå bara i [0, alphaMax].
        inFit = ok;
        t = alphas(inFit)/alphaMax;
        St = alphaMax*(evals.G(:, inFit)'*d);    % lutning dF/dt
        [q, fit, info] = select_poly_degree(t, evals.F(inFit), St);
        evals.fit{end+1} = fit;
        evals.fitInfo{end+1} = info;
        evals.fitAlpha{end+1} = t;
        evals.fitF{end+1} = evals.F(inFit);
        evals.fitSlope{end+1} = St;
        evals.fitAlphaMax{end+1} = alphaMax;
        % Gradvalet under tabellraden (skrivs av print_rest) om showDegree
        degreeLine = '';
        if opts.showDegree
            degreeLine = describe_degree(info);
        end
        evals.Q(end) = q;
        tHat = find_step_length(fit, opts.h, 1);
        alphaHat = tHat*alphaMax;
        evals.AHat(end) = alphaHat;
        if opts.saveFits
            evals.fits{end+1} = struct('k', k, 't', t, 'F', evals.F(inFit), 'St', St, ...
                                       'alphaMax', alphaMax, 'q', q, 'tHat', tHat);
        end

        reason = stop_reason(k, alphas, evals, iBest, alphaHat, alphaMax, d, S0, opts);
        if ~isempty(reason)
            evals.reason = reason;
            print_rest(v, q, [], 'stopp', degreeLine);
            stop = true;
            break
        end
        alphaNew = alphaHat;
        print_rest(v, q, alphaNew, 'polynomets max', degreeLine);
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

    % Ta bort ogiltiga punkter (NaN/Inf) ur evals innan den lämnas
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
    if v
        fprintf('  Stopp efter %d %s%s: %s\n  -> %s\n', n, word, invalidNote, evals.reason, chosen);
    end
end

function reason = stop_reason(k, alphas, evals, iBest, alphaHat, alphaMax, d, S0, opts)
% Stoppkriterierna efter polynomanpassningen. Tomt reason: fortsätt.
% iBest är kolumnen i evals med störst F, alphaHat polynomets maximum.
% (maxEvals och dF/dalpha(0) <= 0 hanteras i line_search_main.)
    reason = '';

    % Liten lutning i bästa punkten (tidigast efter minEvals punkter)
    if k >= opts.minEvals && abs(evals.G(:, iBest)'*d) <= opts.c*abs(S0)
        reason = sprintf('lutningen i bästa punkten alpha = %.3g är liten, |dF/dalpha| <= %.2g*|dF/dalpha(0)|', ...
                         alphas(iBest), opts.c);
        return
    end

    % Maximum ligger närmare än dMin*alphaMax en redan mätt punkt (eller är
    % exakt den). kalmanObj är deterministisk, så en körning där ger inget
    % nytt.
    [dist, i] = min(abs(alphaHat - alphas));     % närmaste mätta punkt
    if dist <= max(opts.dMin*alphaMax, 1e-12*alphaHat)
        reason = sprintf('polynomets max ligger nära alpha = %.3g, redan utvärderad', alphas(i));
    end
end

function s = describe_degree(info)
% Gradvalet som text, en rad per prövat steg: LR mot gränsen och sigma
% (sigmaG i t-skala) före och efter.
    if isempty(info.LR)
        s = sprintf('          gradval  m = %d, qMax = %d, inget test\n', info.m, info.qMax);
        return
    end
    s = sprintf('          gradval  m = %d, qMax = %d\n', info.m, info.qMax);
    for k = 1:numel(info.LR)
        if info.LR(k) > info.critTest(k)
            verdict = sprintf(' > %.2f  höj   ', info.critTest(k));
        else
            verdict = sprintf('<= %.2f  stanna', info.critTest(k));
        end
        s = [s, sprintf('            %d->%d  LR = %7.2f %s   sF %.2e -> %.2e   sG %.2e -> %.2e\n', ...
                        info.qFrom(k), info.qTo(k), info.LR(k), verdict, ...
                        info.sigmaF(1, k), info.sigmaF(2, k), ...
                        info.sigmaG(1, k), info.sigmaG(2, k))]; %#ok<AGROW>
    end
end

function print_rest(v, q, alphaNext, why, degreeLine)
% Avslutar tabellraden: polynomets grad, nästa alpha och varför.
% Tomma q eller alphaNext skrivs som blanka. degreeLine (valfri) skrivs
% på en egen rad under. v = 0: skriv ingenting.
    if ~v
        return
    end
    qStr = '';
    if ~isempty(q)
        qStr = sprintf('%d', q);
    end
    nextStr = '';
    if ~isempty(alphaNext)
        nextStr = sprintf('%.3e', alphaNext);
    end
    fprintf('  %4s  %-9s  %s\n', qStr, nextStr, why);
    if nargin >= 5
        fprintf('%s', degreeLine);
    end
end

function ok = valid_points(evals)
% Punkter som får användas: ändliga F och G.
    ok = isfinite(evals.F) & all(isfinite(evals.G), 1);
end

function opts = default_opts(opts)
% Fyller i standardvärden. Okända fält ger fel, så att ett felstavat eller
% gammalt namn i config-filen inte ignoreras tyst.
    def = struct('h', 1e-3, 'dMin', 1e-2, 'c', 0.05, 'maxEvals', 20, ...
                 'alphaMin', 0.1, 'minEvals', 1, 'showDegree', 0, 'saveFits', 0, ...
                 'shrinkAlphaMax', 1, 'verbose', 0);
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
