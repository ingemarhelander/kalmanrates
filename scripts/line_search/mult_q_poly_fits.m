
function fits = mult_q_poly_fits(alpha, F, S, qMax)
% mult_q_poly_fits - Anpassa ML-polynom av olika gradtal.
%
% In:
%   alpha   Steglängder (m x 1)
%   F       Funktionsvärden (m x 1)
%   S       Riktningsderivator (m x 1)
%   qMax    Högsta polynomgrad som ska testas
%
% Ut:
%   fits    Cell-array med resultat från poly_fit_ml.
%           fits{q} innehåller resultatet för gradtal q.
%
% Anropar poly_fit_ml en gång för varje gradtal q = 1,...,qMax.
% Skapar därefter en gemensam figur med samtliga polynomanpassningar.

    %% Kontrollera indata

    alpha = alpha(:);
    F = F(:);
    S = S(:);

    m = length(alpha);

    if length(F) ~= m || length(S) ~= m
        error('alpha, F and S must have the same number of observations.');
    end

    if m < 1
        error('At least one measurement point is required.');
    end

    if qMax < 1 || qMax ~= floor(qMax)
        error('qMax must be a positive integer.');
    end

    % Gradtalet måste vara strikt mindre än 2*m
    if qMax >= 2*m
        error('qMax must be strictly smaller than 2*m = %d.', 2*m);
    end

    %% Anropa poly_fit_ml för varje gradtal

    fits = cell(qMax, 1);
    rhoF = zeros(qMax, 1);
    rhoG = zeros(qMax, 1);

    % Sortera mätpunkterna efter alpha
    [alphaSorted, sortIdx] = sort(alpha);

    for q = 1:qMax

        fprintf('Fitting polynomial with degree q = %d...\n', q);

        % false innebär att varje anpassning inte skapar en egen figur
        fits{q} = poly_fit_ml(alpha, F, S, q, false);

        % Sortera residualerna enligt alpha
        rF_sorted = fits{q}.rF(sortIdx);
        rG_sorted = fits{q}.rG(sortIdx);

        % Beräkna lag-1-autokorrelation
        rhoF(q) = lag1_autocorr(rF_sorted);
        rhoG(q) = lag1_autocorr(rG_sorted);

        % Spara autokorrelation i fit-resultatet
        fits{q}.rhoF = rhoF(q);
        fits{q}.rhoG = rhoG(q);

        % Undersök linjär trend i residualerna
        [fits{q}.trendCorrF, fits{q}.pTrendF] = ...
            residual_trend_test(alphaSorted, rF_sorted);

        [fits{q}.trendCorrG, fits{q}.pTrendG] = ...
            residual_trend_test(alphaSorted, rG_sorted);

    end

    %% Förbered gemensam plottning

    alphaPlot = linspace(min(alpha), max(alpha), 300)';

    figure('Name', 'Comparison of ML polynomial fits');

    %% Funktionsvärden F

    subplot(2,1,1);

    plot(alpha, F, 'ko', ...
        'MarkerFaceColor', 'b', ...
        'DisplayName', 'Measurements');

    hold on;

    for q = 1:qMax

        y = fits{q}.y;

        % Evaluera polynomet
        Fplot = polyval(flipud(y), alphaPlot);

        plot(alphaPlot, Fplot, ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('q = %d, \\rho_F = %.3f', q, rhoF(q)));

    end

    xlabel('Step length \alpha');
    ylabel('Function value F');
    title('Comparison of polynomial fits for F');
    legend('Location', 'best');
    grid on;
    hold off;


    %% Riktningsderivator S

    subplot(2,1,2);

    plot(alpha, S, 'ko', ...
        'MarkerFaceColor', 'b', ...
        'DisplayName', 'Measurements');

    hold on;

    for q = 1:qMax

        y = fits{q}.y;

        % Koefficienter för derivatpolynomet
        dy = (1:q)' .* y(2:end);

        % Evaluera derivatpolynomet
        Splot = polyval(flipud(dy), alphaPlot);

        plot(alphaPlot, Splot, ...
            'LineWidth', 1.5, ...
            'DisplayName', sprintf('q = %d, \\rho_G = %.3f', q, rhoG(q)));

    end

    xlabel('Step length \alpha');
    ylabel('Derivative S');
    title('Comparison of polynomial fits for S');
    legend('Location', 'best');
    grid on;
    hold off;


    %% Skriv ut resultat

    fprintf('\nML fitting results:\n');

    fprintf(['q    sigmaF    sigmaG    rhoF    rhoG    ' ...
            'pTrendF    pTrendG    iterations\n']);

    for q = 1:qMax

        fprintf('%d    %.5g    %.5g    %.4f    %.4f    %.4g    %.4g    %d\n', ...
            q, ...
            fits{q}.sigmaF, ...
            fits{q}.sigmaG, ...
            fits{q}.rhoF, ...
            fits{q}.rhoG, ...
            fits{q}.pTrendF, ...
            fits{q}.pTrendG, ...
            fits{q}.iterations);

    end
end