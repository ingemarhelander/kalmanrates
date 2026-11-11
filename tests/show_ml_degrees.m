function qs = show_ml_degrees(alpha, F, S, F_true, name)
% show_ml_degrees - Visar gradvalet med ML och lag-1-autokorrelation.
%
% Anpassar ML-polynom av grad q = 2, 3, ... tills residualerna ser ut som
% brus (max rho <= tau), som i select_poly_degree. Ritar en figur med en
% ruta per provad grad: mätpunkterna, polynomet, den sanna funktionen och
% residualerna som streck (rött = punkten ovanför kurvan, blått = nedanför).
%
%   alpha, F, S  mätningar
%   F_true       sann funktion (funktionshandtag), ritas streckad
%   name         namn på fallet, används som rubrik
%   qs           de grader som provades, sista är vald grad
%
% Kräver att scripts/line_search finns på sökvägen.

    % Färger som syns på både ljus och mörk bakgrund
    red    = [1.0 0.3 0.3];
    blue   = [0.3 0.6 1.0];
    orange = [1.0 0.6 0.0];
    green  = [0.2 0.8 0.2];
    gray   = [0.6 0.6 0.6];

    [alpha, order] = sort(alpha(:));
    F = F(:);  F = F(order);
    S = S(:);  S = S(order);
    m = numel(alpha);

    % Anpassa grad för grad, stanna när residualerna ser ut som brus
    fits = {};  rho = [];  tau = [];  qs = [];
    for q = 2:6
        fits{end+1} = poly_fit_ml(alpha, F, S, q, false); %#ok<AGROW>
        rho(end+1, :) = lag1_autocorr([fits{end}.rF, fits{end}.rG]); %#ok<AGROW>
        tau(end+1) = whiteness_threshold(alpha, q, fits{end}); %#ok<AGROW>
        qs(end+1) = q; %#ok<AGROW>
        if max(rho(end, :)) <= tau(end)
            break
        end
    end
    fprintf('%-40s grader: %s\n', name, strjoin(string(qs), ' -> '));

    % En ruta per grad
    n = numel(qs);
    figure('Name', name, 'Position', [50 100 min(450*n, 1800) 450]);
    alphaPlot = linspace(min(alpha), max(alpha), 300)';
    for k = 1:n
        fit = fits{k};
        p = polyval(flipud(fit.y), alpha);
        ok = max(rho(k, :)) <= tau(k);

        subplot(1, n, k);
        hold on;
        for j = 1:m
            if fit.rF(j) >= 0
                color = red;
            else
                color = blue;
            end
            plot([alpha(j) alpha(j)], [p(j) F(j)], 'Color', color, 'LineWidth', 2);
        end
        hTrue = plot(alphaPlot, F_true(alphaPlot), '--', 'Color', gray, 'LineWidth', 1);
        hFit  = plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), '-', 'Color', green, 'LineWidth', 2);
        hData = plot(alpha, F, 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 6);
        grid on;
        hold off;
        xlabel('Steglängd \alpha');
        if k == 1
            ylabel('F');
            legend([hData hFit hTrue], {'Mätningar', 'ML-polynom', 'Sann funktion'}, 'Location', 'best');
        end

        signs = repmat('-', 1, m);
        signs(fit.rF >= 0) = '+';
        if ok
            verdict = 'OK, brus';
        else
            verdict = 'höj graden';
        end
        title({sprintf('q = %d:  %s', qs(k), signs), ...
               sprintf('\\rho_F = %.2f, \\rho_S = %.2f, \\tau = %.2f  =>  %s', ...
                       rho(k, 1), rho(k, 2), tau(k), verdict)});
    end
    sgtitle(name);
end
