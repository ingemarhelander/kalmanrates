% plotOutOfSample.m
% Visualizes the out-of-sample evaluation: runs the Kalman filter from the
% first optimization day to the last out-of-sample day and plots, per day, how well the
% one-step-ahead prediction fits the market (start parameters vs the
% optimized parameters).
%
% Run kalmanOptimizeMain first; it saves data/kalmanOptimizeResult.mat.
% Without that file only the start parameters are plotted.
% The figure is also saved as docs/outOfSample_<model>.png.
%
% Change plotTenors below to choose which contracts to show in the last row.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
plotTenors = {'1Y', '2Y', '5Y', '10Y'}; % Contracts shown in the last row
smoothDays = 10;                         % Rolling mean in the first two rows

config = readConfig(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));
periods = experimentPeriods(config);
if isempty(periods.outOfSampleStart)
  error('kalmanrates:NoOutOfSample', 'Set out_of_sample_start in kalmanOptimizeMain.config to plot the out-of-sample evaluation.');
end
[K, times, ricAll] = loadKalmanModel(periods.modelFile);

% Parameters and periods: from the saved run when it matches the config.
resultFile = dataPath('kalmanOptimizeResult.mat');
if isfile(resultFile)
  R = load(resultFile);
  if ~isfield(R, 'periods') || ~isfield(R.periods, 'optimizeStart')
    warning('kalmanrates:StaleResult', ...
            '%s is from an older version. Run kalmanOptimizeMain again. Plotting the start parameters only.', resultFile);
    R = [];
  elseif ~strcmp(R.modelFile, periods.modelFile)
    warning('kalmanrates:StaleResult', ...
            '%s is from model %s, not %s. Plotting the start parameters only.', ...
            resultFile, R.modelFile, periods.modelFile);
    R = [];
  end
else
  warning('kalmanrates:NoResult', 'No %s. Run kalmanOptimizeMain first. Plotting the start parameters only.', resultFile);
  R = [];
end

xStart = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
if isempty(R)
  paramSets = {xStart, 'Startparametrar'};
else
  periods = R.periods;
  paramSets = {R.xStart, 'Startparametrar' ; R.xOpt, sprintf('Optimerade (%s)', R.savedAt)};
  fprintf('Using %s from %s\n', resultFile, R.savedAt);
end

% Same days as evaluateOutOfSample: first optimization day to last out-of-sample day.
optimizeIdx = selectDays(times, periods.optimizeStart, periods.optimizeEnd);
outOfSampleIdx = selectDays(times, periods.outOfSampleStart, periods.outOfSampleEnd);
idx = optimizeIdx(1):outOfSampleIdx(end);
KS = subsetKalmanModel(K, idx);
days = floor(times(idx));
isOptimize = ismember(idx, optimizeIdx);
isOutOfSample = ismember(idx, outOfSampleIdx);
outOfSampleStartDay = days(find(isOutOfSample, 1));

% Filter each parameter set.
nSets = size(paramSets, 1);
res = cell(nSets, 1);
fprintf('\n%-40s %14s %14s\n', '', 'Optimering', 'Out-of-sample');
for s = 1:nSets
  res{s} = filterDays(KS, paramSets{s,1});
  fMex = kalmanObj(paramSets{s,1}, KS, 0);
  if abs(sum(res{s}.f) - fMex) > 1e-6 * abs(fMex)
    warning('kalmanrates:FilterMismatch', ...
            'The MATLAB filter (%.4f) and mexKalmanLogLikelihood (%.4f) differ.', sum(res{s}.f), fMex);
  end
  perObs = @(m) sum(res{s}.f(m)) / sum(res{s}.nz(m));
  fprintf('%-40s %14.4f %14.4f   log likelihood per observation\n', paramSets{s,2}, perObs(isOptimize), perObs(isOutOfSample));
  fprintf('%-40s %14.3g %14.3g   medel av standardiserat fel^2 (1 = som modellen väntar sig)\n', '', ...
          mean(res{s}.d2(isOptimize)), mean(res{s}.d2(isOutOfSample)));
end

%% Figure
colors = [0.45 0.45 0.45 ; 0 0.4470 0.7410];
fig = figure('Name', 'Out-of-sample', 'Position', [50 50 1200 900]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
title(tl, sprintf('Optimering %s till %s, out-of-sample %s till %s', ...
      datestr(days(find(isOptimize, 1)), 'yyyy-mm-dd'), datestr(days(find(isOptimize, 1, 'last')), 'yyyy-mm-dd'), ...
      datestr(outOfSampleStartDay, 'yyyy-mm-dd'), datestr(max(days(isOutOfSample)), 'yyyy-mm-dd')));

ax1 = nexttile(tl); hold(ax1, 'on');
for s = 1:nSets
  plot(ax1, days, movmean(res{s}.f ./ res{s}.nz, smoothDays), 'Color', colors(s,:), 'LineWidth', 1.2, ...
       'DisplayName', paramSets{s,2});
end
markOutOfSample(ax1, outOfSampleStartDay, days(end));
ylabel(ax1, 'log likelihood / observation');
title(ax1, sprintf('Log likelihood per observation och dag (%d dagars glidande medel), högre är bättre', smoothDays));
legend(ax1, 'Location', 'best');

ax2 = nexttile(tl); hold(ax2, 'on');
for s = 1:nSets
  plot(ax2, days, movmean(res{s}.d2, smoothDays), 'Color', colors(s,:), 'LineWidth', 1.2, ...
       'DisplayName', paramSets{s,2});
end
yline(ax2, 1, 'k:', 'HandleVisibility', 'off');
set(ax2, 'YScale', 'log');
markOutOfSample(ax2, outOfSampleStartDay, days(end));
ylabel(ax2, '\epsilon^T S^{-1} \epsilon / n_z');
title(ax2, 'Standardiserat prediktionsfel i kvadrat (1 = lika stort som modellen väntar sig)');

ax3 = nexttile(tl); hold(ax3, 'on');
best = res{end};
lineColors = lines(numel(plotTenors));
for i = 1:numel(plotTenors)
  j = find(strcmp(ricAll(1:K.nu), ['EUREON' plotTenors{i} '=']), 1);
  if isempty(j)
    continue;
  end
  plot(ax3, days, 100*best.z(:,j), '.', 'Color', lineColors(i,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
  plot(ax3, days, 100*best.pred(:,j), '-', 'Color', lineColors(i,:), 'DisplayName', plotTenors{i});
end
markOutOfSample(ax3, outOfSampleStartDay, days(end));
ylabel(ax3, 'ränta (%)');
title(ax3, sprintf('Marknad (punkter) och prognos en dag framåt (linjer), %s', paramSets{end,2}));
legend(ax3, 'Location', 'best');
linkaxes([ax1 ax2 ax3], 'x');
for a = [ax1 ax2 ax3]
  datetick(a, 'x', 'yyyy-mm', 'keeplimits');
  grid(a, 'on');
end

[~, modelName] = fileparts(periods.modelFile);
if isempty(modelName), modelName = 'kalmanDebug'; end
pngFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', ['outOfSample_' modelName '.png']);
exportgraphics(fig, pngFile, 'Resolution', 150);
fprintf('\nSaved %s\n', pngFile);

%% ------------------------------ Functions -------------------------------

function markOutOfSample(ax, outOfSampleStart, lastDay)
  yl = ylim(ax);
  p = patch(ax, [outOfSampleStart lastDay lastDay outOfSampleStart], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.6], ...
            'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
  uistack(p, 'bottom'); % Lines on top of the shading
  xline(ax, outOfSampleStart, 'r-', 'out-of-sample', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
  ylim(ax, yl);
end
