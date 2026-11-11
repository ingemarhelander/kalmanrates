function fig = plot_fit_residuals(alpha, F, fit, F_true, titleLines)
% plot_fit_residuals - Ritar en ML-anpassning med residualer i en ny figur.
%
% Mätpunkter (orange), ML-polynomet (grönt), sann funktion (grå, streckad)
% och residualerna som streck från kurvan till punkterna: rött = punkten
% ligger ovanför kurvan, blått = nedanför. Färgerna syns på både ljus och
% mörk bakgrund.
%
%   alpha, F    mätningar
%   fit         anpassning från poly_fit_ml (fit.y, fit.rF)
%   F_true      sann funktion (funktionshandtag)
%   titleLines  rubrik, text eller cell med en rad per element

    red    = [1.0 0.3 0.3];
    blue   = [0.3 0.6 1.0];
    orange = [1.0 0.6 0.0];
    green  = [0.2 0.8 0.2];
    gray   = [0.6 0.6 0.6];

    alpha = alpha(:);  F = F(:);
    p = polyval(flipud(fit.y), alpha);       % polynomet i mätpunkterna
    alphaPlot = linspace(min(alpha), max(alpha), 300)';

    fig = figure;
    hold on;
    for j = 1:numel(alpha)
        if fit.rF(j) >= 0
            color = red;
        else
            color = blue;
        end
        plot([alpha(j) alpha(j)], [p(j) F(j)], 'Color', color, 'LineWidth', 2);
    end
    hTrue = plot(alphaPlot, F_true(alphaPlot), '--', 'Color', gray, 'LineWidth', 1);
    hFit  = plot(alphaPlot, polyval(flipud(fit.y), alphaPlot), '-', 'Color', green, 'LineWidth', 2);
    hData = plot(alpha, F, 'o', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 7);
    grid on;
    hold off;

    xlabel('Steglängd \alpha');
    ylabel('F');
    legend([hData hFit hTrue], {'Mätningar', sprintf('ML-polynom, q = %d', numel(fit.y) - 1), ...
            'Sann funktion'}, 'Location', 'best');
    title(titleLines);
end
