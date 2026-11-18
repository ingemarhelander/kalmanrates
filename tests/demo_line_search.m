function demo_line_search()
% demo_line_search - Kör line_search_main på en påhittad målfunktion.
%
% Målfunktionen har två variabler och brus på f och g. tests/stubs/kalmanObj.m
% ersätter Kalman, så varken data eller mex behövs. Två körningar med samma
% slumptal: standardinställningar (dMin = 0.01) och dMin = 0.001, som ger
% fler punkter innan sökningen stoppar.
%
% En figur per varv visar målfunktionen längs d, punkterna hittills,
% ML-polynomet, gradvalet och nästa steglängd. Sista figuren jämför
% körningarna.
%
% Kör från projektroten med:  demo_line_search

    close all;
    testDir = fileparts(mfilename('fullpath'));
    stubDir = fullfile(testDir, 'stubs');
    addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));
    addpath(stubDir);                                 % attrappen före model/kalmanObj
    cleanup = onCleanup(@() rmpath(stubDir));         % tas bort även vid fel

    % Målfunktion f(x) = h(x1) - 2*x2^2 med brus
    h  = @(t) 8*log(1 + 5*t) - 6*t;
    dh = @(t) 40./(1 + 5*t) - 6;
    sF = 0.05;  sG = 0.1;
    K  = @(x) deal(h(x(1)) - 2*x(2)^2 + sF*randn, ...
                   [dh(x(1)); -4*x(2)] + sG*randn(2, 1));

    x0 = [0; 0];
    d  = [2; 0.5];                                    % sökriktning
    fLine = @(a) h(2*a) - 2*(0.5*a).^2;               % f längs d utan brus

    runs = {'standard', struct();
            'dMin = 0.001', struct('dMin', 1e-3)};
    results = cell(1, size(runs, 1));
    for k = 1:size(runs, 1)
        rng(1);
        [f0, g0] = K(x0);
        fprintf('\n%s:\n', runs{k, 1});
        out = struct('name', runs{k, 1}, 'f0', f0);
        [out.alpha, ~, out.fNew, ~, out.evals] = line_search_main(x0, f0, g0, d, K, [], [], runs{k, 2});
        out.S = d'*out.evals.G;                       % lutning längs d i varje punkt
        results{k} = out;
        plotRounds(out, fLine);
    end
    plotComparison(results, fLine);
end

function plotRounds(out, fLine)
% En figur per varv. Varv k: punkterna 0..k (startpunkten och k nya).
    ev = out.evals;
    n = numel(ev.A);
    ok = isfinite(ev.F) & all(isfinite(ev.G), 1);
    aPlot = linspace(0, max(ev.A)*1.1, 300)';

    for k = 2:n
        idx = find(ok(1:k));
        figure('Name', sprintf('%s, varv %d', out.name, k - 1));
        hold on;
        handles = plot(aPlot, fLine(aPlot), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        labels = {'Sann f längs d'};

        % Gör om gradvalet med samma punkter för att visa vad som provades
        loopText = 'ingen anpassning (ingen punkt bättre än startpunkten)';
        if ~isnan(ev.Q(k))
            q = [];  fit = [];  info = [];
            evalc('[q, fit, info] = select_poly_degree(ev.A(idx), ev.F(idx), out.S(idx));');
            loopText = describeDegreeLoop(info, q);
            handles(end+1) = plot(aPlot, poly_value(fit, aPlot), '-', ...
                                  'Color', [0.2 0.8 0.2], 'LineWidth', 2); %#ok<AGROW>
            labels{end+1} = sprintf('ML-polynom, q = %d', q); %#ok<AGROW>
        end
        handles(end+1) = plot(ev.A(1:k), ev.F(1:k), 'o', 'Color', [1 0.6 0], ...
                              'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 8); %#ok<AGROW>
        labels{end+1} = 'Utvärderade punkter'; %#ok<AGROW>
        text(ev.A(1:k), ev.F(1:k), arrayfun(@(i) sprintf('  %d', i), 0:k - 1, 'UniformOutput', false));

        heading = sprintf('%s, varv %d: %d punkter', out.name, k - 1, numel(idx));
        if k < n
            xline(ev.A(k + 1), '-', 'nästa \alpha', 'Color', [0.3 0.6 1], 'LineWidth', 1.5);
            last = sprintf('Nästa \\alpha = %.3f', ev.A(k + 1));
        else
            [~, j] = min(abs(ev.A - out.alpha));
            handles(end+1) = plot(ev.A(j), ev.F(j), 'p', 'Color', [1 0.3 0.3], ...
                                  'MarkerFaceColor', [1 0.3 0.3], 'MarkerSize', 18); %#ok<AGROW>
            labels{end+1} = 'Returnerad punkt'; %#ok<AGROW>
            last = sprintf('Stopp: %s.  Returnerar \\alpha = %.3f, f - f_0 = %+.2f', ...
                           ev.reason, out.alpha, out.fNew - out.f0);
        end
        hold off;
        grid on;
        xlabel('Steglängd \alpha');
        ylabel('f(x_0 + \alpha d)');
        legend(handles, labels, 'Location', 'best');
        title({heading, ['Gradval: ' loopText], last});
    end
end

function s = describeDegreeLoop(info, q)
% Text som beskriver vilka grader som provades och varför loopen stannade.
    parts = {};
    for k = 1:numel(info.LR)
        if info.LR(k) > info.crit
            parts{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f > %.2f, höj', k + 1, k + 2, info.LR(k), info.crit); %#ok<AGROW>
        else
            parts{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f \\leq %.2f, stanna', k + 1, k + 2, info.LR(k), info.crit); %#ok<AGROW>
        end
    end
    if q == info.qMax
        parts{end+1} = sprintf('taket m - 2 = %d', info.qMax);
    end
    s = sprintf('q = %d', q);
    if ~isempty(parts)
        s = [s, '  (', strjoin(parts, ',  '), ')'];
    end
end

function plotComparison(results, fLine)
% Alla utvärderade punkter för båda körningarna i samma figur.
    colors = {[1 0.6 0], [0.3 0.6 1]};
    marks = {'o', 's'};
    aMax = max(cellfun(@(r) max(r.evals.A), results));
    aPlot = linspace(0, aMax*1.1, 300)';

    figure('Name', 'Jämförelse');
    hold on;
    handles = plot(aPlot, fLine(aPlot), '--', 'Color', [0.6 0.6 0.6]);
    labels = {'Sann f längs d'};
    for k = 1:numel(results)
        r = results{k};
        handles(end+1) = plot(r.evals.A, r.evals.F, marks{k}, 'Color', colors{k}, ...
                              'MarkerFaceColor', colors{k}, 'MarkerSize', 8); %#ok<AGROW>
        labels{end+1} = sprintf('%s: %d nya punkter, f - f_0 = %+.2f', ...
                                r.name, numel(r.evals.A) - 1, r.fNew - r.f0); %#ok<AGROW>
    end
    hold off;
    grid on;
    xlabel('Steglängd \alpha');
    ylabel('f(x_0 + \alpha d)');
    legend(handles, labels, 'Location', 'best');
    title('Standard mot dMin = 0.001');
end
