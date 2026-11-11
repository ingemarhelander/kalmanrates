function demo_line_search()
% demo_line_search - Kör hela line_search_main på en påhittad målfunktion.
%
% Målfunktionen har två variabler och brus på f och g. tests/stubs/kalmanObj.m
% ersätter Kalman, så det behövs varken data eller mex. Samma sökriktning
% och samma slumptal körs med gradvalet 'lag1' och 'lr'.
%
% För varje metod ritas en figur per varv: målfunktionen längs sökriktningen,
% punkterna som utvärderats hittills (numrerade i ordning), ML-polynomet med
% vald grad och nästa steglängd. Sista figuren visar varför sökningen
% stoppade och vilken punkt som returneras. En sista figur jämför metoderna.
%
% Kör från projektroten med:  demo_line_search

    close all;
    testDir = fileparts(mfilename('fullpath'));
    stubDir = fullfile(testDir, 'stubs');
    addpath(fullfile(fileparts(testDir), 'scripts', 'line_search'));
    addpath(stubDir);                                 % attrappen före model/kalmanObj
    cleanup = onCleanup(@() rmpath(stubDir));         % tas bort även vid fel

    % Målfunktion f(x) = h(x1) - 2*x2^2 med brus
    h  = @(t) 6*(1 - exp(-3*t)) - 2.5*t.^2;
    dh = @(t) 18*exp(-3*t) - 5*t;
    sF = 0.05;  sG = 0.1;
    K  = @(x) deal(h(x(1)) - 2*x(2)^2 + sF*randn, ...
                   [dh(x(1)); -4*x(2)] + sG*randn(2, 1));

    x0 = [0; 0];
    d  = [2; 0.5];                                    % sökriktning
    fLine = @(a) h(2*a) - 2*(0.5*a).^2;               % f längs d utan brus

    methodNames = {'lag1', 'lr'};
    results = cell(1, 2);
    for k = 1:2
        rng(1);
        [f0, g0] = K(x0);
        fprintf('\n%s:\n', methodNames{k});
        out = struct();
        [out.alpha, ~, out.fNew, ~, out.evals] = line_search_main(x0, f0, g0, d, K, [], [], methodNames{k});
        out.f0 = f0;
        out.S = d'*out.evals.G;                       % lutning längs d i varje punkt
        results{k} = out;
        plotRounds(out, fLine, methodNames{k});
    end
    plotComparison(results, fLine, methodNames);
end

function plotRounds(out, fLine, method)
% En figur per varv. Varv k: k+1 punkter (startpunkten och k nya).
    ev = out.evals;
    n = numel(ev.A);
    ok = isfinite(ev.F) & all(isfinite(ev.G), 1);
    aPlot = linspace(0, max(ev.A)*1.15, 300)';

    for k = 2:n
        idx = find(ok(1:k));
        figure('Name', sprintf('%s, varv %d', method, k - 1));
        hold on;
        hTrue = plot(aPlot, fLine(aPlot), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        handles = hTrue;  labels = {'Sann f längs d'};

        if ~isnan(ev.Q(k))
            fit = [];
            evalc('fit = poly_fit_ml(ev.A(idx), ev.F(idx), out.S(idx), ev.Q(k), false);');
            handles(end+1) = plot(aPlot, polyval(flipud(fit.y), aPlot), '-', ...
                                  'Color', [0.2 0.8 0.2], 'LineWidth', 2); %#ok<AGROW>
            labels{end+1} = sprintf('ML-polynom, q = %d', ev.Q(k)); %#ok<AGROW>
        end
        handles(end+1) = plot(ev.A(1:k), ev.F(1:k), 'o', 'Color', [1 0.6 0], ...
                              'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 8); %#ok<AGROW>
        labels{end+1} = 'Utvärderade punkter'; %#ok<AGROW>
        text(ev.A(1:k), ev.F(1:k), arrayfun(@(i) sprintf('  %d', i), 0:k - 1, 'UniformOutput', false));

        last = k == n;
        if ~last
            xline(ev.A(k + 1), '-', 'nästa \alpha', 'Color', [0.3 0.6 1], 'LineWidth', 1.5);
            heading = sprintf('%s, varv %d: %d punkter, q = %d, nästa \\alpha = %.3f', ...
                              method, k - 1, numel(idx), ev.Q(k), ev.A(k + 1));
            sub = '';
        else
            [~, j] = min(abs(ev.A - out.alpha));
            handles(end+1) = plot(ev.A(j), ev.F(j), 'p', 'Color', [1 0.3 0.3], ...
                                  'MarkerFaceColor', [1 0.3 0.3], 'MarkerSize', 18); %#ok<AGROW>
            labels{end+1} = 'Returnerad punkt'; %#ok<AGROW>
            heading = sprintf('%s, varv %d: %d punkter, q = %d', method, k - 1, numel(idx), ev.Q(k));
            sub = sprintf('Stopp: %s.  Returnerar \\alpha = %.3f, f - f_0 = %+.2f', ...
                          ev.reason, out.alpha, out.fNew - out.f0);
        end
        hold off;
        grid on;
        xlabel('Steglängd \alpha');
        ylabel('f(x_0 + \alpha d)');
        legend(handles, labels, 'Location', 'best');
        title({heading, sub});
    end
end

function plotComparison(results, fLine, methodNames)
% Alla utvärderade punkter för båda metoderna i samma figur.
    colors = {[1 0.6 0], [0.3 0.6 1]};
    marks = {'o', 's'};
    aMax = max(cellfun(@(r) max(r.evals.A), results));
    aPlot = linspace(0, aMax*1.15, 300)';

    figure('Name', 'Jämförelse');
    hold on;
    handles = plot(aPlot, fLine(aPlot), '--', 'Color', [0.6 0.6 0.6]);
    labels = {'Sann f längs d'};
    for k = 1:numel(results)
        r = results{k};
        handles(end+1) = plot(r.evals.A, r.evals.F, marks{k}, 'Color', colors{k}, ...
                              'MarkerFaceColor', colors{k}, 'MarkerSize', 8); %#ok<AGROW>
        labels{end+1} = sprintf('%s: %d nya punkter, returnerar \\alpha = %.3f', ...
                                methodNames{k}, numel(r.evals.A) - 1, r.alpha); %#ok<AGROW>
    end
    hold off;
    grid on;
    xlabel('Steglängd \alpha');
    ylabel('f(x_0 + \alpha d)');
    legend(handles, labels, 'Location', 'best');
    title('Utvärderade punkter med lag-1 och LR');
end
