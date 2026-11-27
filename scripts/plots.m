% plots.m
% Figure 1: objective value and gradient norm after each outer iteration of
%           kalmanOptimizeMain, on a shared x-axis. Iteration 0 is the start point.
%           One curve per start point (numStarts); the best run is drawn thicker.
% Figure 2: market quotes and the model's one-day-ahead prediction for one day
%           at a time (slider), and histograms of the daily error (RMSE) in the
%           optimization period and the out-of-sample period.
% Run kalmanOptimizeMain first; it saves data/kalmanOptimizeResult.mat when
% it is done, and data/kalmanOptimizeProgress.mat after every iteration.
%
% While kalmanOptimizeMain is still running (also in another MATLAB), or
% after it was stopped with Ctrl+C, the progress file is used: every start
% point so far up to its latest iteration, and figure 2 with the best point
% found so far. To choose yourself:
%   plotsUseProgress = true;  plots   % always the progress file
%   plotsUseProgress = false; plots   % always the finished result

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
if exist('kalmanOptimizeSweepPlotData', 'var')
    plot_sweep_results(kalmanOptimizeSweepPlotData);
    return;
end
burnIn = 20; % First filter days left out of the histograms (start-up transient)

progressFile = dataPath('kalmanOptimizeProgress.mat');
if exist('plotsUseProgress', 'var') && ~isempty(plotsUseProgress)
    fromProgress = logical(plotsUseProgress);
else % the progress file when a run is not finished
    fromProgress = false;
    if isfile(progressFile)
        progressState = load(progressFile, 'finished');
        fromProgress = ~progressState.finished;
    end
end

if fromProgress
    [R, runs, best, status] = load_progress(progressFile);
else
    R = load(dataPath('kalmanOptimizeResult.mat'));
    h = R.history;
    if isfield(h, 'runs')
        runs = h.runs;
        best = h.bestRunIndex;
    else % result saved before multistart: one run
        runs = struct('fvals', h.fvals, 'normOfgrads', h.normOfgrads);
        best = 1;
    end
    [runs.done] = deal(true);
    status = sprintf('Klar körning, sparad %s', R.savedAt);
end
fprintf('%s\n', status);

tl = tiledlayout(figure('Name', 'Konvergens'), 2, 1, 'TileSpacing', 'compact');
title(tl, status);
ax1 = nexttile(tl); hold(ax1, 'on');
ax2 = nexttile(tl); hold(ax2, 'on');
set(ax2, 'YScale', 'log');
colors = lines(numel(runs));
for k = 1:numel(runs)
    if isempty(runs(k).fvals)
        continue; % not started yet
    end
    it = 0:numel(runs(k).fvals)-1;
    name = sprintf('start %d', k);
    width = 1;
    if ~runs(k).done
        name = [name ' (pågår)'];
    end
    if k == best
        name = [name ' (bäst)'];
        width = 2.5;
    end
    plot(ax1, it, runs(k).fvals, 'o-', 'Color', colors(k,:), 'LineWidth', width, 'DisplayName', name);
    plot(ax2, it, runs(k).normOfgrads, 'o-', 'Color', colors(k,:), 'LineWidth', width, 'DisplayName', name);
end
ylabel(ax1, 'målfunktion f');
grid(ax1, 'on');
legend(ax1, 'Location', 'best');
ylabel(ax2, '||g||');
xlabel(ax2, 'yttre iteration');
grid(ax2, 'on');
linkaxes([ax1 ax2], 'x');
drawnow; % show the convergence now; the filter below takes a while
fprintf('Ritar marknad mot modell (kör Kalmanfiltret över alla dagar, tar en stund)...\n');

%% Market vs model per day
[K, times, ricAll] = loadKalmanModel(R.modelFile);
optIdx = selectDays(times, R.periods.optimizeStart, R.periods.optimizeEnd);
oosIdx = selectDays(times, R.periods.outOfSampleStart, R.periods.outOfSampleEnd);
idx = optIdx(1):oosIdx(end);
r = filterDays(subsetKalmanModel(K, idx), R.xOpt);

d = struct(); % kalmanOptimizeMain leaves a vector d in the workspace
[d.tenor, order] = sort(cellfun(@ricTenorYears, ricAll(1:K.nu)));
d.z = 100*r.z(:, order);                                 % %
d.pred = 100*r.pred(:, order);
d.rmse = 1e4*sqrt(mean((r.z - r.pred).^2, 2, 'omitnan')); % bp
d.days = floor(times(idx));
d.isOos = ismember(idx, oosIdx)';
isOpt = ismember(idx, optIdx)';
isOpt(1:burnIn) = false;
d.ylim = [min(d.z(burnIn+1:end,:), [], 'all') max(d.z(burnIn+1:end,:), [], 'all')] + [-0.1 0.1]; % same for all days

fig = figure('Name', 'Marknad mot modell');
tl = tiledlayout(fig, 2, 1);
tl.OuterPosition = [0 0.06 1 0.94]; % room for the slider
axC = nexttile(tl);
axH = nexttile(tl); hold(axH, 'on');
edges = linspace(0, prctile(d.rmse(isOpt | d.isOos), 99), 40); % top 1 % left out
histogram(axH, d.rmse(isOpt), edges, 'Normalization', 'probability', 'DisplayName', 'optimering (in-sample)');
histogram(axH, d.rmse(d.isOos), edges, 'Normalization', 'probability', 'DisplayName', 'out-of-sample');
xlabel(axH, 'RMSE per dag (bp)'); ylabel(axH, 'andel dagar');
legend(axH); grid(axH, 'on');

n = numel(d.days);
s = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.1 0.01 0.8 0.03], ...
              'Min', 1, 'Max', n, 'Value', burnIn+1, 'SliderStep', [1 20]/(n-1));
addlistener(s, 'ContinuousValueChange', @(~, ~) drawDay(axC, d, round(s.Value))); % updates while dragging
drawDay(axC, d, burnIn+1);

%% ------------------------------ Functions -------------------------------

function [R, runs, best, status] = load_progress(file)
% The progress of kalmanOptimizeMain as the fields plots needs from
% kalmanOptimizeResult.mat. The best run is the one with the highest
% objective so far (or the one kalmanOptimizeMain chose when finished).
  if ~isfile(file)
    error('kalmanrates:MissingResult', '%s does not exist. Run kalmanOptimizeMain first.', file);
  end
  P = load(file);
  runs = P.runs;
  started = find(arrayfun(@(r) ~isempty(r.fvals), runs));
  if isempty(started)
    error('kalmanrates:MissingResult', 'No start point has finished an iteration yet.');
  end
  if P.finished && ~isempty(P.bestRunIndex)
    best = P.bestRunIndex;
  else
    fNow = arrayfun(@(r) r.f, runs(started));
    fNow(~isfinite(fNow)) = -Inf;
    [~, i] = max(fNow);
    best = started(i);
  end
  R.modelFile = P.modelFile;
  R.periods = P.periods;
  R.xOpt = P.xBase;
  R.xOpt(P.varInd) = runs(best).x;
  nDone = nnz([runs.done]);
  if P.finished
    state = 'klar';
  else
    state = 'pågår eller avbruten';
  end
  status = sprintf(['Progress (%s): %d av %d startpunkter klara, startad %s, ' ...
                    'uppdaterad %s, bäst hittills start %d (f = %.8g)'], state, nDone, ...
                   P.numStarts, P.startedAt, P.updatedAt, best, runs(best).f);
end

function drawDay(ax, d, t)
  have = ~isnan(d.z(t,:));
  plot(ax, d.tenor(have), d.z(t,have), 'o', d.tenor(have), d.pred(t,have), '-');
  if d.isOos(t), tag = 'out-of-sample'; else, tag = 'in-sample'; end
  title(ax, sprintf('%s (%s), RMSE %.1f bp', datestr(d.days(t), 'yyyy-mm-dd'), tag, d.rmse(t)));
  ylim(ax, d.ylim);
  xlabel(ax, 'löptid (år)'); ylabel(ax, 'ränta (%)');
  legend(ax, 'marknad', 'modell, prognos en dag framåt', 'Location', 'southeast');
  grid(ax, 'on');
end

function y = ricTenorYears(ric)
% 'EUREONSW=' -> 1/52, 'EUREON3M=' -> 0.25, 'EUREON2Y=' -> 2, overnight -> 1/365.
  tok = regexp(ric, 'EUREON(\d*|S)([WMY])=', 'tokens', 'once');
  if isempty(tok)
    y = 1/365;
    return;
  end
  n = str2double(tok{1});
  if isnan(n), n = 1; end
  y = n * struct('W', 7/365, 'M', 1/12, 'Y', 1).(tok{2});
end

function plot_sweep_results(data)
numValues = numel(data.histories);
colors = lines(numValues);
fig = figure('Name', 'Kalmanoptimering - parameterjämförelse');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
axF = nexttile(tl);
hold(axF, 'on');
axG = nexttile(tl);
hold(axG, 'on');
set(axG, 'YScale', 'log');

for i = 1:numValues
    h = data.histories{i};
    runIndex = h.bestRunIndex;
    run = h.runs(runIndex);
    iteration = 0:numel(run.fvals)-1;
    plot(axF, iteration, run.fvals, 'o-', 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', data.labels{i});
    plot(axG, iteration, run.normOfgrads, 'o-', 'Color', colors(i,:), ...
        'LineWidth', 2, 'DisplayName', data.labels{i});
end
ylabel(axF, 'målfunktion f');
grid(axF, 'on');
legend(axF, 'Location', 'best', 'Interpreter', 'none');
ylabel(axG, '||g||');
xlabel(axG, 'yttre iteration');
grid(axG, 'on');
linkaxes([axF axG], 'x');
title(tl, 'Konvergens för bästa startpunkt per konfigurationsvärde', ...
    'Interpreter', 'none');

fig = figure('Name', 'Kalmanoptimering - jämförelse av slutresultat');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
axObjective = nexttile(tl);
bar(axObjective, 1:numValues, data.comparison.Objective, 'FaceColor', [0.2 0.45 0.7]);
set(axObjective, 'XTick', 1:numValues, 'XTickLabel', data.labels);
xtickangle(axObjective, 30);
ylabel(axObjective, 'slutligt optimeringsmål');
grid(axObjective, 'on');

axOos = nexttile(tl);
oosPerDay = data.outOfSamplePerDay;
if any(isfinite(oosPerDay))
    bar(axOos, 1:numValues, oosPerDay, 'FaceColor', [0.3 0.65 0.4]);
    set(axOos, 'XTick', 1:numValues, 'XTickLabel', data.labels);
    xtickangle(axOos, 30);
    ylabel(axOos, 'out-of-sample log-likelihood per day');
    grid(axOos, 'on');
else
    axis(axOos, 'off');
    text(axOos, 0.5, 0.5, 'Out-of-sample-period saknas i konfigurationen.', ...
        'HorizontalAlignment', 'center');
end
title(tl, 'Slutresultat per konfigurationsvärde', 'Interpreter', 'none');
end
