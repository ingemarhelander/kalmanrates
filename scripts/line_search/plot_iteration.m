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
%   Rad 4:   gradval: LR-statistik för varje prövat steg q -> q+1 eller q+2 mot
%            sitt gränsvärde (rött streck), för de senaste 6 anpassningarna.
%
% Med backtracking (line_search) finns inga polynom; då ritas bara rad 1.
% Kräver MATLAB R2020b eller senare (tiledlayout i flikar, xline).
%
% In:   history  från kalmanOptimizeMain (efter att stepSizes uppdaterats)
%       evals    utdata från line_search_main
%       fNew     funktionsvärde efter steget
%       gNew     gradient efter steget
%       h        (1e-3) samma som ls_h. Används bara om evals saknar fitTHat
%                (maximum som line_search_main valde, med grid eller newton);
%                då räknas maximum om med grid.
%       label    (valfri) text som läggs till i flikens namn, t.ex. 'denoise',
%                så att två sökningar i samma iteration går att skilja åt

    if nargin < 5 || isempty(h), h = 1e-3; end
    if nargin < 6 || isempty(label)
        suffix = '';
    else
        suffix = [' ' label];
    end
    maxLR = 6;
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
            tHat = chosen_max(evals, k, h);
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
            clim(axF, [1 nFit]);
            cb.Label.String = 'anpassning';
        end

        % ---- Rad 2-3, höger: dF/dalpha mot alpha ----
        nexttile(tl, 10, [2 3]); hold on; grid on
        for k = 1:nFit
            fit  = evals.fit{k};
            aMax = evals.fitAlphaMax{k};
            tHat = chosen_max(evals, k, h);
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
        idx = max(1, nFit - maxLR + 1):nFit;
        for c = 1:numel(idx)
            k    = idx(c);
            fit  = evals.fit{k};
            info = evals.fitInfo{k};
            m    = numel(evals.fitAlpha{k});
            nexttile(tl, 3*nCol + c); hold on; grid on
            if isempty(info.LR)
                text(0.5, 0.5, sprintf('grad %d: för få punkter (m = %d)', ...
                     numel(fit.y) - 1, m), 'HorizontalAlignment', 'center');
                axis([0 1 0 1]);
            else
                % Varje test har egen gräns: crit för q -> q+1, crit2 för q -> q+2
                steps = 1:numel(info.LR);
                bar(steps, info.LR);
                plot(steps, info.critTest, 'r_', 'MarkerSize', 20, 'LineWidth', 2);
                xticks(steps);
                xticklabels(arrayfun(@(s) sprintf('%d->%d', info.qFrom(s), info.qTo(s)), ...
                                     steps, 'UniformOutput', false));
            end
            title(sprintf('anp. %d: \\sigma_F=%.2g, \\sigma_G=%.2g', ...
                          k, fit.sigmaF, fit.sigmaG), 'FontSize', 8);
            if c == 1, ylabel('LR'); end
        end
    end

    tabs.SelectedTab = tab;                   % visa den nya fliken
    drawnow;
    wait_if_paused(fig);
end

function tHat = chosen_max(evals, k, h)
% Polynomets max i anpassning k (i t), som line_search_main valde det.
% Äldre evals utan fitTHat: räkna om med grid.
    if isfield(evals, 'fitTHat') && numel(evals.fitTHat) >= k
        tHat = evals.fitTHat{k};
    else
        tHat = find_step_length(evals.fit{k}, h, 1);
    end
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
