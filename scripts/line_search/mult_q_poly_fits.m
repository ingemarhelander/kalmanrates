function fits = mult_q_poly_fits(alpha, F, S, opts)
% mult_q_poly_fits - Anpassa ML-polynom av olika gradtal.
%
% In:
%   alpha   Steglängder (m x 1)
%   F       Funktionsvärden (m x 1)
%   S       Riktningsderivator (m x 1)
%   opts    (valfri) struct med följande fält:
%             plot     (false)  true: rita figur med alla anpassningar
%             verbose  (false)  true: skriv ut text och rank-varningar
%             qMin     (2)      lägsta gräns för qMax
%             qCap     (10)     högsta gradtal som prövas
%             minSep   (1e-3)   punkter närmare än minSep*(alphaMax-alphaMin)
%                               räknas som samma punkt när qMax bestäms
%
% Ut:
%   fits    Cell-array (nOk x 1) med resultat från poly_fit_ml.
%           fits{q} innehåller resultatet för gradtal q. Bara grader som
%           gick att anpassa ingår, så numel(fits) kan vara mindre än qMax.
%
% Anropar poly_fit_ml en gång för varje gradtal q = 1,...,qMax. Slutar vid
% första grad där punkterna inte bestämmer alla koefficienter (rankbrist),
% eftersom högre grader då också är underbestämda.

    %% Inställningar

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    %% Kontrollera indata

    alpha = alpha(:);
    F = F(:);
    S = S(:);

    m = length(alpha);

    if m < 1
        error('At least one measurement point is required.');
    end
    if length(F) ~= m || length(S) ~= m
        error('alpha, F and S must have the same number of observations.');
    end

    %% Bestäm högsta gradtal

    % Räkna bara punkter som ligger åtskilda. Nästan sammanfallande punkter
    % ger nästan linjärt beroende kolumner och ingen ny information.
    [alphaSorted, sortIdx] = sort(alpha);
    span = alphaSorted(end) - alphaSorted(1);
    if span > 0
        mEff = 1 + sum(diff(alphaSorted) > opts.minSep*span);
    else
        mEff = 1;
    end

    qMax = max(opts.qMin, min(opts.qCap, mEff - 2));

    if qMax < 1 || qMax ~= floor(qMax)
        error('qMax must be a positive integer.');
    end
    % Gradtalet måste vara strikt mindre än 2*m
    if qMax >= 2*m
        error('qMax must be strictly smaller than 2*m = %d.', 2*m);
    end

    %% Anropa poly_fit_ml för varje gradtal

    fits = cell(qMax, 1);
    nOk = 0;

    % Utan verbose tystas MATLABs rank-varning från \ i poly_fit_ml.
    % Fel från rankkontrollen i poly_fit_ml fångas nedan.
    if ~opts.verbose
        ws = warning('off', 'MATLAB:rankDeficientMatrix');
        cleaner = onCleanup(@() warning(ws));
    end

    for q = 1:qMax

        if opts.verbose
            fprintf('Fitting polynomial with degree q = %d...\n', q);
        end

        try
            fit = poly_fit_ml(alpha, F, S, q);
        catch err
            if strcmp(err.identifier, 'poly_fit_ml:rankDeficient')
                if opts.verbose
                    fprintf('  Grad %d är underbestämd, stannar vid grad %d.\n', q, nOk);
                end
                break
            end
            rethrow(err);
        end

        % Sortera residualerna enligt alpha
        rF_sorted = fit.rF(sortIdx);
        rG_sorted = fit.rG(sortIdx);

        % Lag-1-autokorrelation
        fit.rhoF = lag1_autocorr(rF_sorted);
        fit.rhoG = lag1_autocorr(rG_sorted);

        % Linjär trend i residualerna
        [fit.trendCorrF, fit.pTrendF] = residual_trend_test(alphaSorted, rF_sorted);
        [fit.trendCorrG, fit.pTrendG] = residual_trend_test(alphaSorted, rG_sorted);

        fits{q} = fit;
        nOk = q;
    end

    if nOk == 0
        error('mult_q_poly_fits:noFit', ...
              'Ingen polynomgrad gick att anpassa (%d punkter, %d åtskilda).', m, mEff);
    end
    fits = fits(1:nOk);

    %% Figur

    if opts.plot
        plot_fits(alpha, F, S, fits);
    end

    %% Utskrift

    if opts.verbose
        print_table(fits);
    end
end

function plot_fits(alpha, F, S, fits)
% Gemensam figur med alla anpassningar för F och S.

    nOk = numel(fits);
    alphaPlot = linspace(min(alpha), max(alpha), 300)';

    figure('Name', 'Comparison of ML polynomial fits');

    % Funktionsvärden F
    subplot(2,1,1);
    plot(alpha, F, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Measurements');
    hold on;
    for q = 1:nOk
        y = fits{q}.y;
        plot(alphaPlot, polyval(flipud(y), alphaPlot), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('q = %d, \\rho_F = %.3f', q, fits{q}.rhoF));
    end
    xlabel('Step length \alpha');
    ylabel('Function value F');
    title('Comparison of polynomial fits for F');
    legend('Location', 'best');
    grid on;
    hold off;

    % Riktningsderivator S
    subplot(2,1,2);
    plot(alpha, S, 'ko', 'MarkerFaceColor', 'b', 'DisplayName', 'Measurements');
    hold on;
    for q = 1:nOk
        y = fits{q}.y;
        dy = (1:q)' .* y(2:end);
        plot(alphaPlot, polyval(flipud(dy), alphaPlot), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('q = %d, \\rho_G = %.3f', q, fits{q}.rhoG));
    end
    xlabel('Step length \alpha');
    ylabel('Derivative S');
    title('Comparison of polynomial fits for S');
    legend('Location', 'best');
    grid on;
    hold off;
end

function print_table(fits)
% Tabell över anpassningarna.

    fprintf('\nML fitting results:\n');
    fprintf(['q    sigmaF    sigmaG    rhoF    rhoG    ' ...
             'pTrendF    pTrendG    iterations\n']);
    for q = 1:numel(fits)
        fprintf('%d    %.5g    %.5g    %.4f    %.4f    %.4g    %.4g    %d\n', ...
                q, fits{q}.sigmaF, fits{q}.sigmaG, fits{q}.rhoF, fits{q}.rhoG, ...
                fits{q}.pTrendF, fits{q}.pTrendG, fits{q}.iterations);
    end
end

function opts = default_opts(opts)
% Fyller i standardvärden. Okända fält ger fel.

    def = struct('plot', false, 'verbose', false, ...
                 'qMin', 2, 'qCap', 10, 'minSep', 1e-3);

    unknown = setdiff(fieldnames(opts), fieldnames(def));
    if ~isempty(unknown)
        error('mult_q_poly_fits: okänd inställning: %s', strjoin(unknown, ', '));
    end

    names = fieldnames(def);
    for i = 1:numel(names)
        if ~isfield(opts, names{i})
            opts.(names{i}) = def.(names{i});
        end
    end
end