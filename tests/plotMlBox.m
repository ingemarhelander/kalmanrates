function fig = plotMlBox(alpha, F, S, ttl, fig)
% plotMlBox - Visar vad ML-boxen (ml_poly_main) gör med en datamängd.
%
%   Övre vänster  F och ML-polynomet för varje provad grad (vald grad tjock)
%   Övre höger    S och polynomens derivata
%   Nedre vänster residualerna för vald grad, delade med skattad sigma
%   Nedre höger   rho_F och rho_S för varje provad grad, och tröskeln tau
%
% fig (valfri): rita om i en befintlig figur, t.ex. för animering.
% Kräver att scripts/line_search finns på sökvägen.

    if nargin < 5 || isempty(fig) || ~isvalid(fig)
        fig = figure('Name', ttl);
    end
    clf(fig);

    q = [];  fit = [];  info = [];
    evalc('[q, fit, info] = ml_poly_main(alpha, F, S);');   % utan utskrifter

    % select_poly_degree sorterar efter alpha, residualerna är i den ordningen
    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    a = linspace(min(alpha), max(alpha), 200)';
    nq = numel(info.q);
    colors = lines(nq);

    t = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
    title(t, sprintf('%s   (m = %d, vald q = %d)', ttl, numel(alpha), q));

    % F och ML-polynomen
    ax = nexttile(t);  hold(ax, 'on');  grid(ax, 'on');
    for k = 1:nq
        plot(ax, a, poly_value(info.fits{k}, a), 'Color', colors(k, :), ...
             'LineWidth', 1 + 1.5*(info.q(k) == q), 'DisplayName', sprintf('q = %d', info.q(k)));
    end
    plot(ax, alpha, F, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'F');
    xlabel(ax, '\alpha');  ylabel(ax, 'F');  title(ax, 'ML-polynom');
    legend(ax, 'Location', 'best');

    % S och derivatan. Antar att fit.y är koefficienter för 1, alpha, alpha^2, ...
    ax = nexttile(t);  hold(ax, 'on');  grid(ax, 'on');
    for k = 1:nq
        y = info.fits{k}.y;
        dy = (1:numel(y) - 1)'.*y(2:end);
        plot(ax, a, a.^(0:numel(dy) - 1)*dy, 'Color', colors(k, :), ...
             'LineWidth', 1 + 1.5*(info.q(k) == q));
    end
    plot(ax, alpha, S, 'ko', 'MarkerFaceColor', 'k');
    xlabel(ax, '\alpha');  ylabel(ax, 'S = g''d');  title(ax, 'Derivata');

    % Residualer för vald grad, i enheter av skattad sigma
    ax = nexttile(t);  hold(ax, 'on');  grid(ax, 'on');
    stem(ax, alpha, fit.rF/fit.sigmaF, 'filled', 'DisplayName', 'r_F/\sigma_F');
    stem(ax, alpha, fit.rG/fit.sigmaG, 'DisplayName', 'r_S/\sigma_S');
    xlabel(ax, '\alpha');
    title(ax, sprintf('Residualer, q = %d  (\\sigma_F = %.1e, \\sigma_S = %.1e, %d iter)', ...
          q, fit.sigmaF, fit.sigmaG, fit.iterations));
    legend(ax, 'Location', 'best');

    % Autokorrelation och tröskel per provad grad
    ax = nexttile(t);  hold(ax, 'on');  grid(ax, 'on');
    if nq == 1
        % bar tolkar en enda rad som en vektor, lägg till en tom rad
        bar(ax, [info.q, info.q + 1], [info.rho; NaN NaN]);
        xlim(ax, info.q + [-0.5 0.5]);
    else
        bar(ax, info.q, info.rho);
    end
    plot(ax, info.q, info.tau, 'kd-', 'MarkerFaceColor', 'k');
    yline(ax, 0, 'k:');
    ylim(ax, [-1 1]);  xticks(ax, info.q);
    xlabel(ax, 'q');  title(ax, 'Lag-1-autokorrelation (höjs om max \rho > \tau)');
    legend(ax, {'\rho_F', '\rho_S', '\tau'}, 'Location', 'best');

    drawnow;
end
