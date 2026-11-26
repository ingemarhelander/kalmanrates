function plot_iteration(history, evals, fNew, gNew, h, label)
% plot_iteration - Live-figur efter varje linjesökning i kalmanOptimizeMain.
%
% Varje iteration får en egen flik (Iter 1, Iter 2, ...) i samma fönster, så
% du kan bläddra tillbaka i efterhand. De senaste maxTabs (100) flikarna
% sparas. Den nya fliken väljs automatiskt; vill du titta i en äldre flik
% under körningen, tryck först på "Pausa" (annars hoppar fönstret till
% nästa iteration). "Pausa" stoppar körningen efter att fliken ritats,
% tills du trycker på "Fortsätt".
%
% Layout (4 rader x 6 kolumner):
%   Rad 1:   optimeringens förlopp: F, ||g|| (log) och steglängd alpha.
%   Rad 2-3: till vänster F-F0 mot alpha, till höger dF/dalpha mot alpha.
%            Alla ML-polynom från sökningen ligger i samma axlar, färgade
%            efter anpassningsnummer (den senaste tjockast). Triangeln är
%            polynomets valda maximum. Svarta punkter och tangentstreck är
%            de observerade F-värdena och lutningarna i den senaste
%            anpassningen.
%   Rad 4:   gradval för de senaste 6 anpassningarna. Staplarna är poängen
%            från select_poly_degree för varje prövad grad (lägre är bättre).
%            Den valda graden är orange. Med metoden 'threshold' (som saknar
%            poäng) visas i stället vilka grader som klarade alla gränser (1)
%            eller inte (0).
%
% Med backtracking (line_search) finns inga polynom; då ritas bara rad 1.
% Kräver MATLAB R2020b eller senare (tiledlayout i flikar, xline).
%
% In:   history  från kalmanOptimizeMain (efter att stepSizes uppdaterats)
%       evals    utdata från line_search_main
%       fNew     funktionsvärde efter steget
%       gNew     gradient efter steget
%       h        (1e-3) samma som ls_h, används för att hitta polynomets max
%       label    (valfri) text som läggs till i flikens namn, t.ex. 'denoise',
%                så att två sökningar i samma iteration går att skilja åt

    if nargin < 5 || isempty(h), h = 1e-3; end
    if nargin < 6 || isempty(label)
        suffix = '';
    else
        suffix = [' ' label];
    end
    maxFits = 6;                          % antal anpassningar i rad 4
    maxTabs = 100;                        % högsta antal sparade flikar
    nCol = 6;

    % Backtracking-sökningen (line_search) anpassar inga polynom och har
    % varken fit-fälten eller reason; rita då bara förloppet.
    if ~isfield(evals, 'fit') || ~isfield(evals, 'fitAlphaMax')
        evals.fit = {};
    end
    if ~isfield(evals, 'reason')
        evals.reason = '';
    end
    if ~isfield(evals, 'extensionReason') %#ok<NASGU>
        % inget att göra, fältet används inte här
    end

    iter = numel(history.fvals);
    [fig, tabs] = get_figure();
    tab = uitab(tabs, 'Title', sprintf('Iter %d%s', iter, suffix));
    while numel(tabs.Children) > maxTabs      % släng de äldsta flikarna
        delete(tabs.Children(1));
    end
    tl = tiledlayout(tab, 4, nCol, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Iteration %d%s: %s', iter, suffix, evals.reason), 'Interpreter', 'none');

    % ---- Rad 1: optimeringens förlopp ----
    nexttile(tl, 1, [1 2]);
    plot(0:iter, [history.fvals(:); fNew], '-o'); grid on
    title('F'); xlabel('iteration (0 = start)');

    nexttile(tl, 3, [1 2]);
    semilogy(0:iter, [history.normOfgrads(:); norm(gNew)], '-o'); grid on
    title('||g||'); xlabel('iteration (0 = start)');

    nexttile(tl, 5, [1 2]);
    semilogy(1:numel(history.stepSizes), history.stepSizes, '-o'); grid on
    title('steglängd \alpha'); xlabel('iteration');

    nFit = numel(evals.fit);
    if nFit == 0
        nexttile(tl, 7, [3 nCol]);
        text(0.5, 0.5, 'Ingen polynomanpassning gjordes i den här sökningen', ...
             'HorizontalAlignment', 'center');
        axis off
    else
        F0 = evals.F(1);                  % startpunkten är första kolumnen
        % Varje polynom ritas över alla punkter det anpassades till. Alla
        % giltiga punkter ingår, även de bortom alphaMax (t > 1), så kurvan
        % ska inte sluta vid alphaMax.
        tGrid = @(k) linspace(0, max([1; evals.fitAlpha{k}(:)]), 400)';
        cmap = parula(nFit + 1);          % +1 så sista färgen inte är ljusgul

        % ---- Rad 2-3, vänster: F-F0 mot alpha ----
        axF = nexttile(tl, 7, [2 3]); hold on; grid on
        hl = gobjects(nFit, 1);
        labels = cell(nFit, 1);
        for k = 1:nFit
            fit  = evals.fit{k};
            aMax = evals.fitAlphaMax{k};
            tHat = find_step_length(fit, h, 1);
            lw = 1;
            if k == nFit, lw = 2.5; end
            tt = tGrid(k);
            hl(k) = plot(tt*aMax, poly_value(fit, tt) - F0, '-', ...
                         'Color', cmap(k, :), 'LineWidth', lw);
            plot(tHat*aMax, poly_value(fit, tHat) - F0, 'v', ...
                 'Color', cmap(k, :), 'MarkerFaceColor', cmap(k, :), 'MarkerSize', 8);
            labels{k} = sprintf('anpassning %d (grad %d)', k, numel(fit.y) - 1);
        end
        % Observerade punkter och tangenter i den senaste anpassningen
        aMaxL = evals.fitAlphaMax{nFit};
        a = evals.fitAlpha{nFit}(:)*aMaxL;
        F = evals.fitF{nFit}(:) - F0;
        S = evals.fitSlope{nFit}(:);          % dF/dt
        w = 0.03;                             % halv längd på tangentstrecken (i t)
        for i = 1:numel(a)
            plot(a(i) + w*aMaxL*[-1 1], F(i) + S(i)*w*[-1 1], 'k-');
        end
        ho = plot(a, F, 'ko', 'MarkerFaceColor', 'k');
        title('Alla polynom i sökningen (triangel = valt max)');
        xlabel('\alpha'); ylabel('F - F_0');
        if nFit <= 8
            legend([hl; ho], [labels; {'observationer'}], 'Location', 'best');
        else
            colormap(axF, cmap(1:nFit, :));
            cb = colorbar(axF);
            caxis(axF, [1 nFit]);             % caxis fungerar i alla versioner
            cb.Label.String = 'anpassning';
        end

        % ---- Rad 2-3, höger: dF/dalpha mot alpha ----
        nexttile(tl, 10, [2 3]); hold on; grid on
        for k = 1:nFit
            fit  = evals.fit{k};
            aMax = evals.fitAlphaMax{k};
            tHat = find_step_length(fit, h, 1);
            lw = 1;
            if k == nFit, lw = 2.5; end
            tt = tGrid(k);
            plot(tt*aMax, poly_slope(fit, tt)/aMax, '-', ...
                 'Color', cmap(k, :), 'LineWidth', lw);
            plot(tHat*aMax, poly_slope(fit, tHat)/aMax, 'v', ...
                 'Color', cmap(k, :), 'MarkerFaceColor', cmap(k, :), 'MarkerSize', 8);
        end
        plot(a, S/aMaxL, 'ko', 'MarkerFaceColor', 'k');
        yline(0, ':');
        title('Polynomens derivata och observerade lutningar');
        xlabel('\alpha'); ylabel('dF/d\alpha');

        % ---- Rad 4: gradval för de senaste anpassningarna ----
        idx = max(1, nFit - maxFits + 1):nFit;
        for c = 1:numel(idx)
            k = idx(c);
            nexttile(tl, 3*nCol + c);
            plot_degree_choice(evals.fit{k}, evals.fitInfo{k}, k, numel(evals.fitAlpha{k}));
            if c == 1, ylabel('poäng'); end
        end
    end

    tabs.SelectedTab = tab;                   % visa den nya fliken
    drawnow;
    wait_if_paused(fig);
end

function plot_degree_choice(fit, info, k, m)
% Rita gradvalet för en anpassning: poäng per grad, vald grad markerad.
    hold on; grid on
    qSel = numel(fit.y) - 1;                  % vald grad
    ttl = sprintf('anp. %d: \\sigma_F=%.2g, \\sigma_G=%.2g', k, fit.sigmaF, fit.sigmaG);

    if ~isstruct(info) || ~isfield(info, 'degree')
        text(0.5, 0.5, sprintf('grad %d (m = %d): ingen gradinfo', qSel, m), ...
             'HorizontalAlignment', 'center');
        axis([0 1 0 1]);
        title(ttl, 'FontSize', 8);
        return
    end

    degrees = info.degree(:);
    if isfield(info, 'scores') && any(isfinite(info.scores))
        vals = info.scores(:);                % score-metoden: lägre är bättre
        yl = 'poäng';
    else
        vals = double(info.accepted(:));      % threshold-metoden: 1 = godkänd
        yl = 'godkänd';
    end
    vals(~isfinite(vals)) = 0;                % ogiltiga anpassningar ritas som 0

    b = bar(degrees, vals);
    b.FaceColor = 'flat';
    b.CData = repmat([0.6 0.6 0.6], numel(degrees), 1);
    sel = find(degrees == qSel, 1);
    if ~isempty(sel)
        b.CData(sel, :) = [0.85 0.33 0.10];   % vald grad
    end
    xticks(degrees);
    xlabel(sprintf('grad (m = %d)', m), 'FontSize', 8);
    if ~strcmp(yl, 'poäng')
        ylim([0 1.2]);
    end
    title(ttl, 'FontSize', 8);
end

function [fig, tabs] = get_figure()
% Hittar fönstret eller skapar det (med pausknapp) första gången.
    fig = findobj('Type', 'figure', 'Tag', 'kalmanLineSearchLive');
    if ~isempty(fig)
        fig = fig(1);
        tabs = findobj(fig, 'Type', 'uitabgroup');
        if ~isempty(tabs)
            tabs = tabs(1);
            return
        end
        delete(fig);          % gammalt fönster utan flikar: skapa om det
    end
    fig = figure('Tag', 'kalmanLineSearchLive', 'Name', 'Linjesökning (live)', ...
                 'NumberTitle', 'off', 'Units', 'normalized', ...
                 'Position', [0.05 0.05 0.9 0.85]);
    tabs = uitabgroup(fig, 'Units', 'normalized', 'Position', [0 0.05 1 0.95]);
    uicontrol(fig, 'Style', 'togglebutton', 'String', 'Pausa', ...
              'Units', 'pixels', 'Position', [10 10 90 28], ...
              'Tag', 'pauseButton', 'Callback', @toggle_label);
end

function toggle_label(src, ~)
    if src.Value
        src.String = 'Fortsätt';
    else
        src.String = 'Pausa';
    end
end

function wait_if_paused(fig)
% Står kvar här tills knappen släpps (eller fönstret stängs).
    btn = findobj(fig, 'Tag', 'pauseButton');
    while ~isempty(btn) && isgraphics(btn) && btn.Value
        pause(0.1);
    end
end