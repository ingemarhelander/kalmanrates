% plotYieldCurve.m
% Plots the whole EUR OIS curve as the Kalman model sees it, from the first
% optimization day to the last out-of-sample day.
%
% Top:    the model's zero curve (0 to 10 years) for every day as a heatmap.
% Bottom: the curve on a few chosen days: market quotes (points), the model's
%         one-day-ahead prediction of the same contracts (line), the same
%         prediction from the factors only, without the contract-specific
%         states (dotted), and the model's instantaneous forward curve (dashed),
%         where the ECB steps show. Heatmap and forward curve are factors only.
%         A large gap between the dotted and solid lines means the
%         contract-specific states, not the curve factors, carry the fit.
%
% Everything is the one-day-ahead prediction (state before the day's quotes are
% used), so the out-of-sample part is a real forecast.
%
% Run kalmanOptimizeMain first; it saves data/kalmanOptimizeResult.mat.
% Without that file the start parameters are used.
% The figure is also saved as docs/yieldCurve_<model>.png.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
snapshotDates = {};   % e.g. {'2006-06-01', '2006-12-01'}; empty = chosen automatically
maxYears = 10;        % Longest maturity in the heatmap

config = readConfig(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));
periods = experimentPeriods(config);
if isempty(periods.outOfSampleStart)
  error('kalmanrates:NoOutOfSample', 'Set out_of_sample_start in kalmanOptimizeMain.config.');
end
[K, times, ricAll] = loadKalmanModel(periods.modelFile);

% Parameters: optimized ones from the saved run when it matches the model.
paramVec = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
paramName = 'Startparametrar';
resultFile = dataPath('kalmanOptimizeResult.mat');
if isfile(resultFile)
  R = load(resultFile);
  if isfield(R, 'periods') && isfield(R.periods, 'optimizeStart') && strcmp(R.modelFile, periods.modelFile)
    periods = R.periods;
    paramVec = R.xOpt;
    paramName = sprintf('Optimerade (%s)', R.savedAt);
  else
    warning('kalmanrates:StaleResult', '%s does not match %s. Using the start parameters.', resultFile, periods.modelFile);
  end
else
  warning('kalmanrates:NoResult', 'No %s. Run kalmanOptimizeMain first. Using the start parameters.', resultFile);
end

optimizeIdx = selectDays(times, periods.optimizeStart, periods.optimizeEnd);
outOfSampleIdx = selectDays(times, periods.outOfSampleStart, periods.outOfSampleEnd);
idx = optimizeIdx(1):outOfSampleIdx(end);
KS = subsetKalmanModel(K, idx);
days = floor(times(idx));
outOfSampleStartDay = floor(times(outOfSampleIdx(1)));

r = filterStates(KS, paramVec);

% Contract maturities in years, parsed from the RIC (EUREONSW=, EUREON3M=, EUREON2Y=, ...).
ric = ricAll(1:K.nu);
tenorYears = cellfun(@ricTenorYears, ric);

%% Zero curve for every day
nDaysCurve = min(round(maxYears*365), size(r.Ef, 1) - 1);
T = (7:7:nDaysCurve)';               % maturity in days
zeroCurve = nan(numel(T), KS.nt);
for t = 1:KS.nt
  [~, zeroCurve(:,t)] = modelCurves(r, t, T);
end

%% Snapshot days
if isempty(snapshotDates)
  snapDays = [days(find(ismember(idx, optimizeIdx), 1, 'last')), outOfSampleStartDay, ...
              days(round((find(days == outOfSampleStartDay, 1) + KS.nt)/2)), days(end)];
else
  snapDays = cellfun(@(s) datenum(s, 'yyyy-mm-dd'), snapshotDates);
end
snapT = arrayfun(@(d) find(days >= d, 1), snapDays, 'UniformOutput', false);
snapT = unique([snapT{:}], 'stable');

%% Figure
fig = figure('Name', 'Räntekurva', 'Position', [50 50 1300 950]);
tl = tiledlayout(fig, 2, numel(snapT), 'TileSpacing', 'compact');
title(tl, sprintf('Modellens EUR OIS-kurva, prognos en dag framåt, %s', paramName));

axH = nexttile(tl, [1 numel(snapT)]);
imagesc(axH, days, T/365, 100*zeroCurve);
set(axH, 'YDir', 'normal');
colormap(axH, parula);
clim(axH, prctile(100*zeroCurve(:), [1 99])); % the first days are start-up noise
cb = colorbar(axH); cb.Label.String = 'nollkupongränta (%)';
hold(axH, 'on');
xline(axH, outOfSampleStartDay, 'r-', 'out-of-sample', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
for k = snapT
  xline(axH, days(k), 'w--', 'HandleVisibility', 'off');
end
datetick(axH, 'x', 'yyyy-mm', 'keeplimits');
ylabel(axH, 'löptid (år)');
title(axH, 'Nollkupongkurva per dag (vita streck = dagarna nedan)');

for k = snapT
  ax = nexttile(tl); hold(ax, 'on');
  [fwd, ~] = modelCurves(r, k, (1:nDaysCurve)');
  plot(ax, (1:nDaysCurve)/365, 100*fwd, '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'modell, forward');
  have = ~isnan(r.z(k,:));
  [tY, order] = sort(tenorYears(have));
  zk = r.z(k, have); pk = r.pred(k, have);
  sk = r.syst(k, have);
  plot(ax, tY, 100*sk(order), ':', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.4, ...
       'DisplayName', 'modell, bara faktorer');
  plot(ax, tY, 100*pk(order), '-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.4, 'DisplayName', 'modell, swapränta');
  plot(ax, tY, 100*zk(order), 'o', 'Color', [0.85 0.325 0.098], 'MarkerFaceColor', [0.85 0.325 0.098], ...
       'MarkerSize', 4, 'DisplayName', 'marknad');
  rmse = 1e4*sqrt(mean((zk - pk).^2));
  if days(k) >= outOfSampleStartDay, tag = 'out-of-sample'; else, tag = 'optimering'; end
  title(ax, sprintf('%s (%s)\nRMSE %.1f bp', datestr(days(k), 'yyyy-mm-dd'), tag, rmse));
  xlabel(ax, 'löptid (år)');
  xlim(ax, [0 maxYears]);
  grid(ax, 'on');
  if k == snapT(1)
    ylabel(ax, 'ränta (%)');
    legend(ax, 'Location', 'best');
  end
end

[~, modelName] = fileparts(periods.modelFile);
pngFile = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'docs', ['yieldCurve_' modelName '.png']);
exportgraphics(fig, pngFile, 'Resolution', 150);
fprintf('Saved %s\n', pngFile);

%% ------------------------------ Functions -------------------------------

function r = filterStates(KS, paramVec)
% Kalman filter as in plotOutOfSample (clsKalmanIr.logLikelihoodEps). Stores per
% day the predicted systematic state, the observed quotes z and the
% one-day-ahead prediction pred of each contract.

  n = cumsum([length(KS.thetaF) length(KS.Sigmav) length(KS.Sigmaw) length(KS.ax) length(KS.Sigma0)]);
  thetaF = paramVec(1:n(1));
  Sigmav = paramVec(n(1)+1:n(2));
  Sigmaw = paramVec(n(2)+1:n(3));
  ax0    = paramVec(n(3)+1:n(4));
  Sigma0 = paramVec(n(4)+1:n(5));
  Ef = KS.Ef; Ef(:) = paramVec(n(5)+1:end);
  EfInt = [zeros(1, size(Ef,2)) ; cumsum(Ef,1)]/365;

  nt = KS.nt;
  r.Ef = Ef; r.EfInt = EfInt; r.np = KS.np; r.nc = KS.nc;
  r.datesStepRel = KS.datesStepRelAll;
  r.xs = zeros(KS.ns, nt);
  r.z = nan(nt, KS.nu); r.pred = nan(nt, KS.nu); r.syst = nan(nt, KS.nu);
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
    datesStepRel = KS.datesStepRelAll{t};
    o = [-EfInt(KS.oTAll{t},:) zeros(size(KS.oTAll{t},1), KS.nc)];
    for i = 1:KS.nc
      o(:,KS.np+i) = -max(KS.oTAll{t}-(datesStepRel(i)+1), 0)/365;
    end
    [g, G] = taylorInstrPC(KS.pl, KS.atAll{t}, o, KS.oIndAll{t}, KS.tcAll{t}, xPred(1:KS.ns));
    u = g - G*xPred(1:KS.ns);
    H = [G KS.Iz{t}];
    zPred = H*xPred + u;
    S = H*PPred*H' + Cv;
    Kg = (PPred*H')/S;
    xUpdate = xPred + Kg*(KS.z{t} - zPred);
    PUpdate = PPred - Kg*H*PPred;

    r.xs(:,t) = xPred(1:KS.ns);
    [row, col] = find(KS.Iz{t});
    r.z(t, col) = KS.z{t}(row);
    r.pred(t, col) = zPred(row);
    r.syst(t, col) = g(row);   % factors only, without the contract-specific states
  end
end

function [fwd, zero] = modelCurves(r, t, T)
% Instantaneous forward rate and zero rate at maturities T (days) on day t:
% f(T) = Ef(T)*xp + sum_i xc_i*[T >= step_i], zero(T) = integral of f / T.
  xp = r.xs(1:r.np, t);
  xc = r.xs(r.np+1:end, t);
  dsr = r.datesStepRel{t}(:)';
  fwd  = r.Ef(T+1,:)*xp    + double(T >= dsr)*xc;
  zero = (r.EfInt(T+1,:)*xp + max(T - dsr, 0)/365*xc) ./ (T/365);
end

function y = ricTenorYears(ric)
% 'EUREONSW=' -> 1/52, 'EUREON3M=' -> 0.25, 'EUREON2Y=' -> 2.
  tok = regexp(ric, 'EUREON(\d*|S)([WMY])=', 'tokens', 'once');
  if isempty(tok)
    y = NaN;
    return;
  end
  n = str2double(tok{1});
  if isnan(n), n = 1; end
  y = n * struct('W', 7/365, 'M', 1/12, 'Y', 1).(tok{2});
end
