% plotOutOfSample.m
% Visualizes the out-of-sample test: runs the Kalman filter from the first
% training day to the last test day and plots, per day, how well the
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
if isempty(periods.testStart)
  error('kalmanrates:NoTest', 'Set test_start in kalmanOptimizeMain.config to plot the out-of-sample test.');
end
[K, times, ricAll, KI] = loadKalmanModel(periods.modelFile);

% Parameters, days and contracts: from the saved run when it matches the config.
resultFile = dataPath('kalmanOptimizeResult.mat');
paramSets = {};
if isfile(resultFile)
  R = load(resultFile);
  if ~strcmp(R.modelFile, periods.modelFile)
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
  trainOpts = struct('startDate', periods.trainStart, 'endDate', periods.trainEnd, ...
      'weekdays', cfgValue(config, 'weekdays', 'all'), 'dayStep', cfgValue(config, 'day_step', 1));
  testOpts = struct('startDate', periods.testStart, 'endDate', periods.testEnd);
  keepInstr = selectInstruments(ricAll, K.nu, cfgValue(config, 'tenors', 'all'));
  paramSets = {xStart, 'Startparametrar'};
else
  trainOpts = R.dataSel.train;
  testOpts = R.dataSel.test;
  keepInstr = R.keepInstr;
  paramSets = {R.xStart, 'Startparametrar' ; R.xOpt, sprintf('Optimerade (%s)', R.savedAt)};
  fprintf('Using %s from %s\n', resultFile, R.savedAt);
end

% Same days as evaluateOutOfSample: first training day to last test day.
trainIdx = selectDays(times, trainOpts);
filters = trainOpts;
filters.startDate = floor(times(trainIdx(1)));
filters.endDate = testOpts.endDate;
idx = selectDays(times, filters);
KS = subsetKalmanModel(K, idx, keepInstr);
days = floor(times(idx));
isTrain = ismember(idx, trainIdx);
isTest = days >= testOpts.startDate & days <= testOpts.endDate;
testStartDay = min(days(isTest));

% Filter each parameter set.
nSets = size(paramSets, 1);
res = cell(nSets, 1);
fprintf('\n%-40s %14s %14s\n', '', 'Träning', 'Test');
for s = 1:nSets
  res{s} = filterDays(KS, paramSets{s,1});
  fMex = kalmanObj(paramSets{s,1}, KS, 0);
  if abs(sum(res{s}.f) - fMex) > 1e-6 * abs(fMex)
    warning('kalmanrates:FilterMismatch', ...
            'The MATLAB filter (%.4f) and mexKalmanLogLikelihood (%.4f) differ.', sum(res{s}.f), fMex);
  end
  perObs = @(m) sum(res{s}.f(m)) / sum(res{s}.nz(m));
  fprintf('%-40s %14.4f %14.4f   log likelihood per observation\n', paramSets{s,2}, perObs(isTrain), perObs(isTest));
  fprintf('%-40s %14.3g %14.3g   medel av standardiserat fel^2 (1 = som modellen väntar sig)\n', '', ...
          mean(res{s}.d2(isTrain)), mean(res{s}.d2(isTest)));
end

%% Figure
colors = [0.45 0.45 0.45 ; 0 0.4470 0.7410];
fig = figure('Name', 'Out-of-sample', 'Position', [50 50 1200 900]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
title(tl, sprintf('Träning %s till %s, test %s till %s', ...
      datestr(days(find(isTrain, 1)), 'yyyy-mm-dd'), datestr(days(find(isTrain, 1, 'last')), 'yyyy-mm-dd'), ...
      datestr(testStartDay, 'yyyy-mm-dd'), datestr(max(days(isTest)), 'yyyy-mm-dd')));

ax1 = nexttile(tl); hold(ax1, 'on');
for s = 1:nSets
  plot(ax1, days, movmean(res{s}.f ./ res{s}.nz, smoothDays), 'Color', colors(s,:), 'LineWidth', 1.2, ...
       'DisplayName', paramSets{s,2});
end
markTest(ax1, testStartDay, days(end));
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
markTest(ax2, testStartDay, days(end));
ylabel(ax2, '\epsilon^T S^{-1} \epsilon / n_z');
title(ax2, 'Standardiserat prediktionsfel i kvadrat (1 = lika stort som modellen väntar sig)');

ax3 = nexttile(tl); hold(ax3, 'on');
[~, labels] = selectInstruments(ricAll, K.nu, 'all');
best = res{end};
lineColors = lines(numel(plotTenors));
for i = 1:numel(plotTenors)
  j = find(strcmp(labels, plotTenors{i}) & keepInstr, 1);
  if isempty(j)
    continue;
  end
  plot(ax3, days, 100*best.z(:,j), '.', 'Color', lineColors(i,:), 'MarkerSize', 6, 'HandleVisibility', 'off');
  plot(ax3, days, 100*best.pred(:,j), '-', 'Color', lineColors(i,:), 'DisplayName', plotTenors{i});
end
markTest(ax3, testStartDay, days(end));
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

function r = filterDays(KS, paramVec)
% Kalman filter as clsKalmanIr.logLikelihoodEps, for the struct that
% subsetKalmanModel returns. Per day: log likelihood f, number of
% observations nz, standardized squared error d2 = eps'*inv(S)*eps/nz, and
% per contract the observed z and the one-step-ahead prediction pred.

  thetaFInd = (1:length(KS.thetaF));
  SigmavInd = (1:length(KS.Sigmav)) + thetaFInd(end);
  SigmawInd = (1:length(KS.Sigmaw)) + SigmavInd(end);
  axInd     = (1:length(KS.ax))     + SigmawInd(end);
  Sigma0Ind = (1:length(KS.Sigma0)) + axInd(end);
  EfVecInd  = (1:length(KS.Ef(:)))  + Sigma0Ind(end);
  thetaF = paramVec(thetaFInd);
  Sigmav = paramVec(SigmavInd);
  Sigmaw = paramVec(SigmawInd);
  ax0    = paramVec(axInd);
  Sigma0 = paramVec(Sigma0Ind);
  Ef = KS.Ef; Ef(:) = paramVec(EfVecInd);
  EfInt = [zeros(1, size(Ef,2)) ; cumsum(Ef,1)]/365;

  nt = KS.nt;
  r.f = zeros(nt,1); r.nz = zeros(nt,1); r.d2 = zeros(nt,1);
  r.z = nan(nt, KS.nu); r.pred = nan(nt, KS.nu);
  for t = 1:nt
    F  = sparse(KS.A{t}*diag(thetaF)*KS.B{t});
    Cv = KS.G{t}*diag(Sigmav)*KS.G{t}';
    Cw = KS.D{t}*diag(Sigmaw)*KS.D{t}';
    if t == 1
      xPred = F*ax0;
      PPred = F*diag(Sigma0)*F' + Cw;
    else
      xPred = F*xUpdate;
      PPred = F*PUpdate*F' + Cw;
    end
    o = oCompute(KS, t, EfInt);
    [g, G] = taylorInstrPC(KS.pl, KS.atAll{t}, o, KS.oIndAll{t}, KS.tcAll{t}, xPred(1:KS.ns));
    u = g - G*xPred(1:KS.ns);
    H = [G KS.Iz{t}];
    zPred = H*xPred + u;
    eps = KS.z{t} - zPred;
    S = H*PPred*H' + Cv;
    Kg = (PPred*H')/S;
    xUpdate = xPred + Kg*eps;
    PUpdate = PPred - Kg*H*PPred;

    w = S \ eps;
    r.nz(t) = KS.nz(t);
    r.f(t) = -0.5*(KS.nz(t)*log(2*pi) + log(det(S)) + eps'*w);
    r.d2(t) = (eps'*w) / KS.nz(t);
    [row, col] = find(KS.Iz{t});
    r.z(t, col) = KS.z{t}(row);
    r.pred(t, col) = zPred(row);
  end
end

function o = oCompute(KS, t, EfInt)
% Same as the private clsKalmanIr.oCompute.
  datesStepRel = KS.datesStepRelAll{t};
  o = [-EfInt(KS.oTAll{t},:) zeros(size(KS.oTAll{t},1), KS.nc)];
  for i = 1:KS.nc
    o(:,KS.np+i) = -max(KS.oTAll{t}-(datesStepRel(i)+1), 0)/365;
  end
end

function markTest(ax, testStart, lastDay)
  yl = ylim(ax);
  p = patch(ax, [testStart lastDay lastDay testStart], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.6], ...
            'FaceAlpha', 0.35, 'EdgeColor', 'none', 'HandleVisibility', 'off');
  uistack(p, 'bottom'); % Lines on top of the shading
  xline(ax, testStart, 'r-', 'test', 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
  ylim(ax, yl);
end

function value = cfgValue(config, key, default)
  value = default;
  if isfield(config, key)
    value = config.(key);
  end
end
