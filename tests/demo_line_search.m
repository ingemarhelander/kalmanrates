function demo_line_search(method)
% demo_line_search - Kör hela line_search_main på en påhittad målfunktion.
%
% Målfunktionen har två variabler och brus på f och g. tests/stubs/kalmanObj.m
% ersätter Kalman, så det behövs varken data eller mex. Line search körs två
% gånger med samma slumptal:
%   1. med standardinställningarna,
%   2. med minPoints = 8: sökningen får inte stoppa på lutning eller brus
%      förrän 8 punkter finns, så att gradvalet hinner prova fler grader.
%
% För varje körning ritas en figur per varv: målfunktionen längs
% sökriktningen, punkterna hittills (numrerade i ordning), ML-polynomet och
% nästa steglängd. Titeln visar gradvalet i varvet: vilka grader som
% provades med de punkter som finns och varför det stannade. Sista figuren
% i varje körning visar stoppskälet och den returnerade punkten. En sista
% figur jämför körningarna.
%
% method (valfri): 'lr' (standard) eller 'lag1'.
% Kör från projektroten med:  demo_line_search

    if nargin < 1
        method = 'lr';
    end
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
            'minPoints = 8', struct('minPoints', 8)};
    results = cell(1, size(runs, 1));
    for k = 1:size(runs, 1)
        rng(1);
        [f0, g0] = K(x0);
        fprintf('\n%s, %s:\n', method, runs{k, 1});
        out = struct('name', runs{k, 1}, 'method', method, 'f0', f0);
        opts = runs{k, 2};
        opts.degreeMethod = method;
        [out.alpha, ~, out.fNew, ~, out.evals] = line_search_main(x0, f0, g0, d, K, [], [], opts);
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
    qPrev = [];

    for k = 2:n
        idx = find(ok(1:k));
        figure('Name', sprintf('%s, varv %d', out.name, k - 1));
        hold on;
        handles = plot(aPlot, fLine(aPlot), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        labels = {'Sann f längs d'};

        % Gör om gradvalet med samma punkter för att visa vad som provades
        loopText = 'ingen anpassning (för få giltiga punkter)';
        if ~isnan(ev.Q(k))
            q = [];  fit = [];  info = [];
            evalc('[q, fit, info] = ml_poly_main(ev.A(idx), ev.F(idx), out.S(idx), qPrev, out.method);');
            loopText = describeDegreeLoop(out.method, info, q, numel(idx));
            qPrev = q;
            handles(end+1) = plot(aPlot, polyval(flipud(fit.y), aPlot), '-', ...
                                  'Color', [0.2 0.8 0.2], 'LineWidth', 2); %#ok<AGROW>
            labels{end+1} = sprintf('ML-polynom, q = %d', q); %#ok<AGROW>
        end
        handles(end+1) = plot(ev.A(1:k), ev.F(1:k), 'o', 'Color', [1 0.6 0], ...
                              'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 8); %#ok<AGROW>
        labels{end+1} = 'Utvärderade punkter'; %#ok<AGROW>
        text(ev.A(1:k), ev.F(1:k), arrayfun(@(i) sprintf('  %d', i), 0:k - 1, 'UniformOutput', false));

        heading = sprintf('%s (%s), varv %d: %d punkter', out.method, out.name, k - 1, numel(idx));
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

function s = describeDegreeLoop(method, info, q, m)
% Text som beskriver vilka grader som provades och varför loopen stannade.
    parts = {};
    if strcmp(method, 'lr')
        for k = 1:numel(info.LR)
            if info.LR(k) > info.crit
                parts{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f > %.2f, höj', k + 1, k + 2, info.LR(k), info.crit); %#ok<AGROW>
            else
                parts{end+1} = sprintf('%d\\rightarrow%d: LR = %.1f \\leq %.2f, stanna', k + 1, k + 2, info.LR(k), info.crit); %#ok<AGROW>
            end
        end
        if q == info.qMax
            parts{end+1} = sprintf('taket m - 2 = %d', info.qMax); %#ok<AGROW>
        end
    else
        qMax = 2;
        if m >= 4
            qMax = min(6, 2*m - 3);
        end
        for k = 1:numel(info.q)
            if m < 4
                parts{end+1} = 'färre än 4 punkter: bara q = 2'; %#ok<AGROW>
            elseif max(info.rho(k, :)) <= info.tau(k)
                parts{end+1} = sprintf('q = %d: \\rho = %.2f \\leq \\tau = %.2f, OK', info.q(k), max(info.rho(k, :)), info.tau(k)); %#ok<AGROW>
            else
                parts{end+1} = sprintf('q = %d: \\rho = %.2f > \\tau = %.2f, höj', info.q(k), max(info.rho(k, :)), info.tau(k)); %#ok<AGROW>
            end
        end
        if m >= 4 && q == qMax
            parts{end+1} = sprintf('taket %d', qMax); %#ok<AGROW>
        end
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
    title(sprintf('%s: standard mot minPoints = 8', results{1}.method));
end
