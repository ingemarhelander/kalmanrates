function plot_polynomial_extension(evals, extendedEvals, ...
                                   extendedFit, d, iter)
% Visar ursprungligt och utökat ML-polynom i en egen flik.
% Använder samma figur och fliksystem som plot_iteration.
% Påverkar inte optimeringen.

    if isempty(evals.fit) || isempty(extendedFit)
        return
    end

    maxTabs = 100;

    % Hämta samma figur och flikgrupp som plot_iteration använder.
    [fig, tabs] = get_extension_figure();

    % Skapa en ny flik i samma fönster.
    tab = uitab(tabs, ...
        'Title', sprintf('Iter %d extension', iter));

    % Behåll högst maxTabs flikar.
    while numel(tabs.Children) > maxTabs
        delete(tabs.Children(1));
    end

    % Senaste ursprungliga och utökade fit.
    originalFit = evals.fit{end};
    originalScale = evals.fitAlphaMax{end};
    extendedScale = extendedEvals.fitAlphaMax{end};

    F0 = evals.F(1);

    % Rita över intervallet från noll till största utvärderade alpha.
    alphaMaxPlot = max([1.3; extendedEvals.A(:)]);
    alphaPlot = linspace(0, alphaMaxPlot, 500)';

    tOriginal = alphaPlot / originalScale;
    tExtended = alphaPlot / extendedScale;

    fOriginal = poly_value(originalFit, tOriginal) - F0;
    fExtended = poly_value(extendedFit, tExtended) - F0;

    % Identifiera vilka punkter som är nya.
    isNew = true(size(extendedEvals.A));

    for k = 1:numel(extendedEvals.A)
        isNew(k) = ~any(abs(evals.A - extendedEvals.A(k)) < 1e-10);
    end

    % Layout
    tl = tiledlayout(tab, 2, 2, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    title(tl, sprintf('Iteration %d: experimentell polynomutvidgning', iter));

    % ---- Vänster: funktionsvärde ----
    ax1 = nexttile(tl, 1, [2 1]);
    hold(ax1, 'on');
    grid(ax1, 'on');

    plot(ax1, alphaPlot, fOriginal, '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Ursprungligt fit');

    plot(ax1, alphaPlot, fExtended, '-', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Utökat fit');

    scatter(ax1, evals.A, evals.F - F0, 45, 'b', 'filled', ...
        'DisplayName', 'Ursprungliga punkter');

    scatter(ax1, extendedEvals.A(isNew), ...
        extendedEvals.F(isNew) - F0, 55, 'g', 'filled', ...
        'DisplayName', 'Nya punkter');

    xline(ax1, 0.95, ':', 'alpha = 0.95');
    xline(ax1, 1.0, ':', 'alpha = 1.0');

    xlabel(ax1, '\alpha');
    ylabel(ax1, 'F - F_0');
    title(ax1, 'Funktionsvärde och polynomfit');
    legend(ax1, 'Location', 'best');

    % ---- Höger upptill: derivata ----
    ax2 = nexttile(tl, 2);
    hold(ax2, 'on');
    grid(ax2, 'on');

    slopeOriginal = poly_slope(originalFit, tOriginal) / originalScale;
    slopeExtended = poly_slope(extendedFit, tExtended) / extendedScale;

    plot(ax2, alphaPlot, slopeOriginal, '--', ...
        'LineWidth', 1.5, ...
        'DisplayName', 'Ursprungligt fit');

    plot(ax2, alphaPlot, slopeExtended, '-', ...
        'LineWidth', 2.5, ...
        'DisplayName', 'Utökat fit');

    yline(ax2, 0, ':');

    xlabel(ax2, '\alpha');
    ylabel(ax2, 'dF/d\alpha');
    title(ax2, 'Polynomets derivata');
    legend(ax2, 'Location', 'best');

    % ---- Höger nedtill: information ----
    ax3 = nexttile(tl, 4);
    axis(ax3, 'off');

    originalDegree = numel(originalFit.y) - 1;
    extendedDegree = numel(extendedFit.y) - 1;

    if isfield(extendedEvals, 'extensionReason')
        reason = extendedEvals.extensionReason;
    else
        reason = 'Okänd';
    end

    if isfield(extendedEvals, 'extensionEvaluations')
        nNew = extendedEvals.extensionEvaluations;
    else
        nNew = nnz(isNew);
    end

    infoText = {
        sprintf('Ursprunglig grad: %d', originalDegree)
        sprintf('Utökad grad: %d', extendedDegree)
        sprintf('Nya Kalman-utvärderingar: %d', nNew)
        sprintf('Totalt antal punkter: %d', numel(extendedEvals.A))
        ''
        'Stopporsak:'
        reason
    };

    text(ax3, 0.05, 0.95, infoText, ...
        'Units', 'normalized', ...
        'VerticalAlignment', 'top', ...
        'Interpreter', 'none');

    % Visa den nya fliken.
    tabs.SelectedTab = tab;
    drawnow;

    wait_if_extension_paused(fig);
end


function [fig, tabs] = get_extension_figure()
% Använder samma figur-tag som plot_iteration.

    fig = findobj('Type', 'figure', ...
                  'Tag', 'kalmanLineSearchLive');

    if ~isempty(fig)
        fig = fig(1);
        tabs = findobj(fig, 'Type', 'uitabgroup');

        if ~isempty(tabs)
            tabs = tabs(1);
            return
        end

        delete(fig);
    end

    % Om extension anropas innan plot_iteration skapat fönstret.
    fig = figure( ...
        'Tag', 'kalmanLineSearchLive', ...
        'Name', 'Linjesökning (live)', ...
        'NumberTitle', 'off', ...
        'Units', 'normalized', ...
        'Position', [0.05 0.05 0.9 0.85]);

    tabs = uitabgroup(fig, ...
        'Units', 'normalized', ...
        'Position', [0 0.05 1 0.95]);

    uicontrol(fig, ...
        'Style', 'togglebutton', ...
        'String', 'Pausa', ...
        'Units', 'pixels', ...
        'Position', [10 10 90 28], ...
        'Tag', 'pauseButton', ...
        'Callback', @toggle_extension_label);
end


function toggle_extension_label(src, ~)
    if src.Value
        src.String = 'Fortsätt';
    else
        src.String = 'Pausa';
    end
end


function wait_if_extension_paused(fig)
    btn = findobj(fig, 'Tag', 'pauseButton');

    while ~isempty(btn) && isgraphics(btn) && btn.Value
        pause(0.1);
    end
end
