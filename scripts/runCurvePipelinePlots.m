% Smooth the states with the optimized parameters and plot the
% interest-rate curve.
% Similar setup to runReducedIrCurvesV2.m, but with additional plotting of the
%smoothed states and reconstructed curves.
%
% Nothing is built or optimized again. The flow is
%   buildKalmanModel; kalmanOptimizeMain; runCurvePipelinePlots
% This script loads the model buildKalmanModel saved for the periods in
% kalmanOptimizeMain.config and the parameters kalmanOptimizeMain saved in
% data/kalmanOptimizeResult.mat.
% The state plots mark the optimization and out-of-sample periods.
% Set rebuildMarketCurves to false to reuse the existing EUR100.mat curves.


%-------------- Load Optimized Parameters ---------------------------
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
run(fullfile(repoDir, 'setupPaths.m'));

% false: use the existing EUR100.mat, the same curves buildKalmanModel used.
% true rebuilds every day of the history in parfor, which needs a lot of
% memory (mexPortfolio "bad allocation" on a laptop) and can make the model
% differ from the one kalmanOptimizeMain optimized.
rebuildMarketCurves = false;
showCurveMovie = false;
showStepComponentMovie = true;

runStart = tic; % Elapsed time in the progress prints: [   12.3 s]

if rebuildMarketCurves
  fprintf('[%7.1f s] Rebuilding the market curves (runTenorMultiIrCurves)...\n', toc(runStart));
  run(fullfile(scriptDir, 'runTenorMultiIrCurves.m')); %This script builds the EUR100.mat curves and saves them to dataPath('EUR100.mat').
end

% Periods from kalmanOptimizeMain.config.
periods = experimentPeriods(readConfig(fullfile(scriptDir, 'kalmanOptimizeMain.config')));
if periods.hasInsample
  % The model buildKalmanModel saved for these periods: K (clsKalmanIr with
  % the start parameters), KI and P (curves and ECB dates for the plots).
  modelPath = dataPath(periods.modelFile);
  if ~isfile(modelPath)
    error('kalmanrates:MissingModel', 'The model %s does not exist. Run buildKalmanModel first.', modelPath);
  end
  fprintf('[%7.1f s] Loading the model %s...\n', toc(runStart), modelPath);
  saved = load(modelPath, 'K', 'KI', 'P');
  if ~isfield(saved, 'P')
    error('kalmanrates:MissingModel', ...
          '%s was built before the plot data was saved. Run buildKalmanModel again.', modelPath);
  end
  K = saved.K;
  KI = saved.KI;
  cb = saved.P.cb;
  cbEffectiveDates = saved.P.cbEffectiveDates;
  nSteps = saved.P.nSteps;
  flgTick = saved.P.flgTick;
  fH = saved.P.fH;
  times = KI.times;
  indInstrAll = KI.indInstrAll;
  clear saved
else
  % Without in-sample there is no saved model with plot data, so build it
  % as runReducedIrCurvesV2 does.
  fprintf('[%7.1f s] No in-sample in the config: building the model (runReducedIrCurvesInit/Cls)...\n', toc(runStart));
  run(fullfile(scriptDir, 'runReducedIrCurvesInit.m'));
  run(fullfile(scriptDir, 'runReducedIrCurvesCls.m'));
end
fprintf('[%7.1f s] Model ready: %d days, %s to %s, %d states\n', toc(runStart), K.nt, ...
        datestr(floor(times(1)), 'yyyy-mm-dd'), datestr(floor(times(end)), 'yyyy-mm-dd'), K.nx);

% Optimized parameters from the last kalmanOptimizeMain run.
resultFile = dataPath('kalmanOptimizeResult.mat');
if ~isfile(resultFile)
  error('kalmanrates:NoResult', 'No %s. Run kalmanOptimizeMain first.', resultFile);
end
optimizerResult = load(resultFile, 'xOpt', 'modelFile', 'savedAt');
if ~strcmp(optimizerResult.modelFile, periods.modelFile)
  error('kalmanrates:StaleResult', ...
        '%s is from model "%s" but the config gives "%s". Run kalmanOptimizeMain with the current config.', ...
        resultFile, optimizerResult.modelFile, periods.modelFile);
end
optimizedParameters = optimizerResult.xOpt;
expectedParameterCount = numel(K.thetaF) + numel(K.Sigmav) + ...
    numel(K.Sigmaw) + numel(K.ax) + numel(K.Sigma0) + numel(K.Ef);
if numel(optimizedParameters) ~= expectedParameterCount
  error('kalmanrates:ModelAlignment', ...
        'The optimized parameters have %d values but the rebuilt model needs %d.', ...
        numel(optimizedParameters), expectedParameterCount);
end
fprintf('[%7.1f s] Using the optimized parameters from %s (%s)\n', toc(runStart), resultFile, optimizerResult.savedAt);
if periods.hasInsample
  [~, runOutputSuffix] = fileparts(periods.modelFile);
  runOutputSuffix = ['_' runOutputSuffix];
end
K.storeParamVec(optimizedParameters);
paramVec = optimizedParameters;

% ---- RTS state-smoothin from runReducedCurvesV2.m------
nRtsIter = 5;
fprintf('[%7.1f s] RTS smoothing (at most %d iterations; columns: f before, f after, change)\n', toc(runStart), nRtsIter);
for iter = 1:nRtsIter
  fprintf('[%7.1f s]   RTS iteration %d/%d...\n', toc(runStart), iter, nRtsIter);
  [KR, f0, fInit] = kalmanSub(K);
  if iter <= 100
    lambda = 1;
  elseif iter < 10
    lambda = 0.01;
  else
    lambda = 1;
  end
  K.x0 = K.x0 + lambda * KR.dx0;
  for t = 1:K.nt
    K.x{t} = K.x{t} + lambda * KR.dx{t};
    K.w{t} = K.w{t} + lambda * KR.dw{t};
    K.v{t} = K.v{t} + lambda * KR.dv{t};
  end
  fNew = K.evalObj();
  fprintf('            %.16d %.16d %d\n', f0, fNew, fNew - f0);
  if abs(f0 - fInit) < 1E-8
    fprintf('[%7.1f s]   RTS converged after %d iterations\n', toc(runStart), iter);
    break;
  end
end

Ef = K.Ef;
axRtsOpt = K.x{1};
axRtsFile = dataPath(['axRtsOpt' runOutputSuffix '.mat']);
save(axRtsFile, 'axRtsOpt');
fprintf('[%7.1f s] Saved %s\n', toc(runStart), axRtsFile);

% Diagnostic methods defined in model/clsKalmanIr.m.
fprintf('[%7.1f s] K.plotSolution (figures 1-6)...\n', toc(runStart));
K.plotSolution(KI, paramVec);
fprintf('[%7.1f s] K.analyze (figures 11-16)...\n', toc(runStart));
K.analyze(KI, paramVec);

% State extraction and component plots copied from
% runReducedIrCurvesV2.m, lines 44-88.
xp = zeros(K.nt, K.np);
xc = zeros(K.nt, K.nc);
xu = zeros(K.nt, K.nu);
v = zeros(K.nt, K.nu);
w = zeros(K.nt, K.nx);
z = nan(K.nt, K.nx);
for t = 1:K.nt
  xp(t,:) = K.x{t}(1:K.np);
  xc(t,:) = K.x{t}(K.np + (1:K.nc));
  xu(t,:) = K.x{t}(K.ns+1:end);
  v(t,indInstrAll{t}) = K.v{t};
  w(t,:) = K.w{t};
  z(t,indInstrAll{t}) = K.z{t};
end

%% Figure
% Against dates, with the optimization and out-of-sample periods marked.
fprintf('[%7.1f s] State plots (figures 21-24)...\n', toc(runStart));
stateFigures = {21, xp(:,1:min(3,K.np)), 'x PC' ; 22, xc(:,1:3), 'x step' ; ...
                23, xu, 'x unsystematic' ; 24, v, 'v'};
for i = 1:size(stateFigures, 1)
  figure(stateFigures{i,1}); clf;
  plot(floor(times), stateFigures{i,2});
  title(stateFigures{i,3});
  datetick('x', 'yyyy-mm', 'keeplimits');
  grid on;
  markPeriods(gca, periods, floor(times(end)));
end

% ECB policy-rate plot and PDF export copied from
% runReducedIrCurvesV2.m, lines 91-111.
fprintf('[%7.1f s] ECB policy-rate plot (figure 25)...\n', toc(runStart));
ind = find(~isnan(cb{:,7}));
fInc = zeros(size(times));
cbRate = zeros(size(times));
for i = 2:length(ind)
  fInc(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = ...
      cb{ind(i),7} - cb{ind(i-1),7};
  cbRate(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = ...
      cb{ind(i),7};
end
ecbFig = figure(25);
clf(ecbFig);
% Whole model period when the config chooses the days, otherwise from 2022.
ind = times >= datenum(2022,01,01);
if periods.hasInsample || ~any(ind)
  ind = true(size(times));
end
plot(times(ind), cbRate(ind) + xc(ind,1)*100, times(ind), cbRate(ind));
datetick('x', 'yyyy');
title('x step');
legend('Market implied policy rate', 'ECB policy rate', 'AutoUpdate', 'off')
markPeriods(gca, periods, floor(times(end)));
% exportgraphics as in plotOutOfSample: crops to the content and also works
% for docked figures, where saveas can fail with "Invalid graphics object".
ecbPdf = fullfile(repoDir, 'docs', ['ecbCutsHikes' runOutputSuffix '.pdf']);
exportgraphics(ecbFig, ecbPdf, 'ContentType', 'vector');
fprintf('[%7.1f s] Saved %s\n', toc(runStart), ecbPdf);

pause(0.01);

% Curve movie copied from runReducedIrCurvesV2.m, lines 117-149.
if showCurveMovie % So the user can see the evolution of the curves over time.
  fprintf('[%7.1f s] Curve movie (figure 26%s), %d days...\n', toc(runStart), ...
          repmat(' and 7-8', 1, showStepComponentMovie), K.nt);
  figure(26);
  for k = 1:K.nt
    if mod(k, 50) == 0
      fprintf('[%7.1f s]   day %d/%d (%s)\n', toc(runStart), k, K.nt, datestr(floor(times(k)), 'yyyy-mm-dd'));
    end
    tradeDate = floor(times(k));
    datesStep = cbEffectiveDates(cbEffectiveDates > tradeDate);
    datesStep = datesStep(1:nSteps);

    Es = zeros(size(Ef,1), nSteps);
    for i = 1:nSteps
      Es(datesStep(i)-tradeDate+1:end, i) = 1;
    end
    E = [Ef Es];

    if (flgTick)
      plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]');
      title(datestr(times(k)));
    else
      plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]', ...
          (0:size(fH,2)-1)/365, fH(k,:));
      title(datestr(floor(times(k))));
    end
    pause(0.01);

    % Optional step-component movie from the if (false) block in
    % runReducedIrCurvesV2.m; use separate figures 7-8.
    if showStepComponentMovie
      figure(7 + mod(k,2));
      plot((0:size(E,1)-1)/365, E(:,K.np+1:end)*[xc(k,:)]');
      title(datestr(floor(times(k))));
      hold off;

      figure(7 + mod(k+1,2));
      plot((0:size(E,1)-1)/365, E(:,K.np+1:end)*[xc(k,:)]');
      hold on;
    end
  end
else
  fprintf('[%7.1f s] Curve movie skipped (showCurveMovie = false)\n', toc(runStart));
end
fprintf('[%7.1f s] runCurvePipelinePlots done\n', toc(runStart));

%% ------------------------------ Functions -------------------------------

function markPeriods(ax, periods, lastDay)
% Shades the optimization period (blue) and the out-of-sample period
% (orange) from kalmanOptimizeMain.config. Nothing is drawn for a period
% that is not set.
  hold(ax, 'on');
  yl = ylim(ax);
  spans = {periods.optimizeStart, periods.optimizeEnd, [0.6 0.75 1], 'optimering' ; ...
           periods.outOfSampleStart, min([periods.outOfSampleEnd lastDay]), [1 0.85 0.6], 'out-of-sample'};
  for i = 1:size(spans, 1)
    [first, last, color, name] = spans{i,:};
    if isempty(first) || isempty(last)
      continue;
    end
    p = patch(ax, [first last last first], [yl(1) yl(1) yl(2) yl(2)], color, ...
              'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    uistack(p, 'bottom'); % Lines on top of the shading
    xline(ax, first, '-', name, 'Color', color*0.6, 'HandleVisibility', 'off', ...
          'LabelVerticalAlignment', 'bottom');
  end
  ylim(ax, yl);
  hold(ax, 'off');
end
