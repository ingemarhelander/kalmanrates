% cleanMarketData.m
% Removes quotes that deviate on their own from the market data, so that the
% model can be built and optimized again on cleaned data (pass 2).
%
% Needs pass 1: buildKalmanModel and kalmanOptimizeMain on the raw data
% (data/kalmanOptimizeResult.mat and the model file it refers to).
%
% 1. Filters the days of the model, from the first optimization day, with
%    the optimized parameters (innovationsLOO). A quote is flagged when its
%    leave-one-out innovation, its error given the past and all other quotes
%    the same day, exceeds the threshold. Rate changes and other days when
%    many contracts move together are explained by the other quotes and are
%    not flagged.
% 2. Prints checks: share of flagged quotes, flags on ECB meeting days, days
%    when many contracts moved (should keep their quotes), contracts that are
%    often flagged, and plots some flagged quotes.
% 3. Writes data/EUR_clean.mat: EUR_raw.mat with BID = NaN for the flagged
%    quotes, and makes it the active data set (useDataSet('clean')). The
%    original is kept as data/EUR_raw.mat and data/EUR100_raw.mat.
%    The flags are saved in data/noisyQuotes.mat.
%
% Why only BID: populatePortfolioWithData drops days where both BID and ASK
% are NaN, and mexPortfolio then returns the quote of the day before (one
% day old, so it is still used). With only BID NaN the day is stored and the
% quote is skipped that day by runTenorMultiIrCurves and buildKalmanModel.
%
% Then, pass 2 (no existing file is changed):
%   runTenorMultiIrCurves   % optional: curves (EUR100.mat) of the cleaned data,
%                           % so that the PCA is built on cleaned curves. Slow.
%   buildKalmanModel
%   kalmanOptimizeMain
%   comparePasses           % checks the cleaning and compares pass 1 and 2
%
% useDataSet('raw') / useDataSet('clean') switches between the data sets.
% To only look at what would be flagged:  dryRun = true; cleanMarketData

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));

%% Settings
% dryRun = true: only flag, print the checks and plot. Nothing is written
% and the active data set is unchanged. Can be set before the script is run:
%   dryRun = true; cleanMarketData
if ~exist('dryRun', 'var') || isempty(dryRun)
  dryRun = false;
end
threshold = 4;            % flag when |rLOO| > threshold (rLOO is N(0,1) when the model holds)
maxFlagFraction = 0.2;    % remove at most this share of the quotes of one day
keepPreviousFlags = true; % keep the flags of an earlier run (noisyQuotes.mat), e.g. after pass 2
moveFraction = 0.3;       % "market move day": at least this share of the quotes has |rStd| > threshold
calibrationRange = [0.5 2]; % robust std of rLOO must lie here, otherwise the parameters do not fit the data
allowUncalibrated = false;  % true: write the cleaned data even when the calibration check fails
nExamplePlots = 6;        % flagged quotes to plot
plotWindow = 15;          % days before and after a flagged quote in the plots
resultFile = 'kalmanOptimizeResult.mat';

%% 1. Pass 1: parameters and model
R = load(dataPath(resultFile));
fprintf('Pass 1: %s (saved %s), model %s\n', resultFile, R.savedAt, R.modelFile);
[K, times, ricAll, KI] = loadKalmanModel(R.modelFile);
if isempty(R.periods.optimizeStart)
  firstDay = 1;
else
  firstDay = selectDays(times, R.periods.optimizeStart, []);
  firstDay = firstDay(1);
end
days = firstDay:K.nt;
KS = subsetKalmanModel(K, days);
dayTimes = floor(times(days));
fprintf('Filtering %d days, %s to %s, %d contracts\n', numel(days), ...
        datestr(dayTimes(1), 'yyyy-mm-dd'), datestr(dayTimes(end), 'yyyy-mm-dd'), K.nu);

%% 2. Raw data and contracts that cannot be removed
% The raw data is EUR_raw.mat once it exists, before that the active
% EUR.mat (which then is the raw data). It is saved as EUR_raw.mat only
% when the cleaned data is written (step 6).
rawFile = dataPath('EUR_raw.mat');
if ~isfile(rawFile)
  if ~strcmp(useDataSet(), 'raw')
    error('kalmanrates:MissingDataSet', 'EUR_raw.mat is missing although the active data set is not raw.');
  end
  rawFile = dataPath('EUR.mat');
end
raw = load(rawFile);
% Only BID/ASK quotes can be removed. The overnight fixing (IBOR sheet) is
% always used by the model and is never flagged.
protected = ~ismember(ricAll, raw.ric{1});
if any(protected)
  fprintf('Never removed (not a BID/ASK quote): %s\n', strjoin(ricAll(protected), ', '));
end

%% 3. Flag quotes
res = innovationsLOO(KS, R.xOpt, struct('threshold', threshold, ...
                     'maxFlagFraction', maxFlagFraction, 'protected', protected));

[ti, ci] = find(res.flag);
ind = sub2ind(size(res.flag), ti, ci);
flags = table(dayTimes(ti), reshape(ricAll(ci), [], 1), res.z(ind), res.pred(ind), ...
              res.rRemoved(ind), res.rLOO(ind), res.rStd(ind), ...
              'VariableNames', {'date', 'ric', 'z', 'pred', 'rRemoved', 'rLOO', 'rStd'});
flags.source = repmat({R.modelFile}, height(flags), 1);
nNew = height(flags);

%% 4. Checks
quoted = ~isnan(res.z);
nObs = nnz(quoted);

% Calibration: rLOO is N(0,1) only when the parameters fit the data. The
% robust std (1.4826*MAD) measures the bulk without the outliers. Far from 1
% the threshold means nothing: much smaller and nothing can be flagged
% (the model believes the noise is larger than it is), much larger and
% ordinary quotes are flagged.
robustStd = @(v) 1.4826*mad(v(~isnan(v)), 1);
bulkStd = robustStd(res.rLOO(quoted));
calibrated = bulkStd >= calibrationRange(1) && bulkStd <= calibrationRange(2);
fprintf('\n=== Calibration ===\n');
fprintf('Robust std of rLOO: %.3g (about 1 when the parameters fit the data)\n', bulkStd);
contractStd = arrayfun(@(i) robustStd(res.rLOO(:,i)), 1:K.nu);
off = find(contractStd < calibrationRange(1) | contractStd > calibrationRange(2));
if ~isempty(off) && calibrated
  fprintf('Contracts outside %.1f-%.1f (their flags are less reliable):\n', calibrationRange);
  fprintf('  %s\n', strjoin(arrayfun(@(i) sprintf('%s %.2g', ricAll{i}, contractStd(i)), off, ...
                                     'UniformOutput', false), ', '));
end
if ~calibrated
  fprintf(['NOT CALIBRATED: the parameters of pass 1 do not describe the size of the\n' ...
           'innovations (%.3g instead of about 1), so the flags below cannot be trusted.\n' ...
           'Optimize until Sigmav and Sigmaw fit the data, then run cleanMarketData again.\n'], bulkStd);
  svInd = numel(K.thetaF) + (1:numel(K.Sigmav));
  fprintf('  median Sigmav %.3g (variance of the quote noise); daily change of the quotes, std %.3g\n', ...
          median(R.xOpt(svInd)), std(reshape(diff(res.z), [], 1), 'omitnan'));
end
fprintf('\n=== Flagged quotes (|rLOO| > %g) ===\n', threshold);
fprintf('%d of %d quotes (%.3f %%) on %d of %d days\n', nNew, nObs, 100*nNew/nObs, ...
        nnz(res.nFlag), numel(days));
expected = nObs * erfc(threshold/sqrt(2));
fprintf('Expected by chance if the model held: about %.1f\n', expected);
if nNew/nObs > 0.01
  fprintf('NOTE: more than 1 %% flagged. The model may not describe the data well;\n');
  fprintf('      look at the plots before using the cleaned data.\n');
end
if any(res.capped)
  fprintf('%d days reached maxFlagFraction (%.0f %%), check them by hand:\n', nnz(res.capped), 100*maxFlagFraction);
  fprintf('  %s\n', strjoin(cellstr(datestr(dayTimes(res.capped), 'yyyy-mm-dd')), ', '));
end

% ECB days: meeting day, the day after, and the day the new rate takes effect
% (as in runReducedIrCurvesInit)
try
  cb = readtable(dataPath('raw', 'datesEURECB.xlsx'), 'Sheet', 'mHistory', 'Format', 'auto');
  meetings = floor(datenum(cb.MeetingDates(~isnat(cb.MeetingDates))));
  effective = datenum(cb.StartDateImplied_assumeThursday_);
  actual = datenum(cb.StartDate_MainRefinancing);
  effective(~isnan(actual)) = actual(~isnan(actual));
  effective = floor(effective(~isnan(effective)));
  ecbDay = ismember(dayTimes, meetings) | ismember(dayTimes - 1, meetings) | ismember(dayTimes, effective);
  rateEcb = nnz(res.flag(ecbDay,:)) / max(nnz(quoted(ecbDay,:)), 1);
  rateOther = nnz(res.flag(~ecbDay,:)) / max(nnz(quoted(~ecbDay,:)), 1);
  fprintf('\nECB days (meeting, day after, effective date): %d days\n', nnz(ecbDay));
  fprintf('  flagged share on ECB days %.3f %%, other days %.3f %%\n', 100*rateEcb, 100*rateOther);
  ecbFlagDays = dayTimes(ecbDay & res.nFlag > 0);
  if ~isempty(ecbFlagDays)
    fprintf('  ECB days with flags: %s%s\n', strjoin(cellstr(datestr(ecbFlagDays(1:min(10, end)), 'yyyy-mm-dd')), ', '), ...
            repmat(' ...', 1, numel(ecbFlagDays) > 10));
  end
  if rateEcb > 2*rateOther && nnz(res.flag(ecbDay,:)) >= 3
    fprintf('  NOTE: clearly more flags on ECB days. Real rate moves may be flagged;\n');
    fprintf('        look at those days before using the cleaned data.\n');
  end
catch err
  fprintf('\nECB dates could not be read (%s), the ECB check is skipped.\n', err.message);
end

% Days when many contracts moved: the ordinary innovation is large for many
% quotes. These should keep their quotes.
bigStd = abs(res.rStd) > threshold;
moveDay = sum(bigStd, 2) >= moveFraction * res.nz;
fprintf('\nMarket move days (at least %.0f %% of the quotes with |rStd| > %g): %d\n', ...
        100*moveFraction, threshold, nnz(moveDay));
if any(moveDay)
  fprintf('  quotes with |rStd| > %g on those days: %d, of which flagged: %d\n', threshold, ...
          nnz(bigStd(moveDay,:)), nnz(res.flag(moveDay,:)));
  list = find(moveDay);
  list = list(1:min(10, end));
  fprintf('  %s%s\n', strjoin(cellstr(datestr(dayTimes(list), 'yyyy-mm-dd')), ', '), ...
          repmat(' ...', 1, nnz(moveDay) > 10));
end

% Contracts that are often flagged (possibly illiquid)
flagRate = sum(res.flag, 1) ./ max(sum(quoted, 1), 1);
[~, order] = sort(flagRate, 'descend');
fprintf('\nContracts with the highest share of flagged quotes:\n');
for i = order(1:min(5, end))
  if flagRate(i) == 0, break; end
  fprintf('  %-14s %5.2f %% (%d of %d)\n', ricAll{i}, 100*flagRate(i), nnz(res.flag(:,i)), nnz(quoted(:,i)));
end
if any(flagRate > 0.05)
  fprintf('  NOTE: a contract with more than 5 %% flagged may be better left out entirely.\n');
end

if nNew > 0
  [~, order] = sort(abs(flags.rRemoved), 'descend');
  fprintf('\nLargest deviations:\n');
  fprintf('  %-10s %-14s %10s %10s %9s %9s\n', 'date', 'ric', 'z', 'pred', 'rRemoved', 'rStd');
  for i = reshape(order(1:min(10, end)), 1, [])
    fprintf('  %-10s %-14s %10.5f %10.5f %9.2f %9.2f\n', datestr(flags.date(i), 'yyyy-mm-dd'), ...
            flags.ric{i}, flags.z(i), flags.pred(i), flags.rRemoved(i), flags.rStd(i));
  end
end

%% Plots
figure('Name', 'cleanMarketData: leave-one-out innovations');
v = res.rLOO(quoted);
edges = linspace(-max(8, threshold+2), max(8, threshold+2), 81);
histogram(v, edges, 'Normalization', 'pdf');
hold on;
xx = linspace(edges(1), edges(end), 400);
plot(xx, exp(-xx.^2/2)/sqrt(2*pi), 'LineWidth', 1.5);
xline([-threshold threshold], '--');
set(gca, 'YScale', 'log');
legend('rLOO', 'N(0,1)', 'threshold');
title(sprintf('Leave-one-out innovations, %d quotes, %d flagged', nObs, nNew));
xlabel('rLOO');

if nNew > 0
  figure('Name', 'cleanMarketData: flagged quotes');
  [~, order] = sort(abs(flags.rRemoved), 'descend');
  nPlot = min(nExamplePlots, nNew);
  for p = 1:nPlot
    i = order(p);
    t = find(dayTimes == flags.date(i), 1);
    c = find(strcmp(ricAll, flags.ric{i}), 1);
    w = max(1, t - plotWindow):min(numel(days), t + plotWindow);
    subplot(ceil(nPlot/2), 2, p);
    plot(dayTimes(w), res.z(w, c), '.-', dayTimes(w), res.pred(w, c), '-');
    hold on;
    plot(dayTimes(t), res.z(t, c), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    datetick('x', 'yyyy-mm-dd', 'keeplimits');
    title(sprintf('%s %s, rLOO %.1f', flags.ric{i}, datestr(flags.date(i), 'yyyy-mm-dd'), flags.rRemoved(i)));
    if p == 1, legend('quote', 'prediction', 'flagged', 'Location', 'best'); end
  end
end

if dryRun
  fprintf('\ndryRun: nothing is written. Set dryRun = false (or clear dryRun) to write the cleaned data.\n');
  return;
end
if ~calibrated && ~allowUncalibrated
  error('kalmanrates:Uncalibrated', ...
        ['The robust std of rLOO is %.3g, outside %.1f-%.1f, so nothing is written. Optimize pass 1 ' ...
         'further, or set allowUncalibrated = true to write anyway.'], ...
        bulkStd, calibrationRange(1), calibrationRange(2));
end

% Keep the raw data before anything is written
if strcmp(useDataSet(), 'raw')
  useDataSet('backup');
end

%% 5. Combine with earlier flags
flagFile = dataPath('noisyQuotes.mat');
if keepPreviousFlags && isfile(flagFile)
  prev = load(flagFile, 'flags');
  key = @(T) strcat(cellstr(datestr(T.date, 'yyyy-mm-dd')), '|', T.ric);
  old = prev.flags(~ismember(key(prev.flags), key(flags)), :);
  fprintf('\nKeeping %d flags from the earlier run (%s)\n', height(old), flagFile);
  flags = [old ; flags];
end
flags = sortrows(flags, {'date', 'ric'});

%% 6. Write the cleaned data
clean = raw;
nMarked = 0;
missing = {};
for i = 1:height(flags)
  row = find(clean.dates{1} == flags.date(i), 1);
  rowAsk = find(clean.dates{2} == flags.date(i), 1);
  col = find(strcmp(clean.ric{1}, flags.ric{i}), 1);
  colAsk = find(strcmp(clean.ric{2}, flags.ric{i}), 1);
  if isempty(row) || isempty(rowAsk) || isempty(col) || isempty(colAsk)
    missing{end+1} = sprintf('%s %s', datestr(flags.date(i), 'yyyy-mm-dd'), flags.ric{i}); %#ok<SAGROW>
    continue;
  end
  % BID = NaN removes the quote. ASK must not be NaN, otherwise the day is
  % dropped and the quote of the day before is used (0 is only a marker; a
  % quote with BID = NaN is never used).
  clean.data{1}(row, col) = NaN;
  if isnan(clean.data{2}(rowAsk, colAsk))
    clean.data{2}(rowAsk, colAsk) = 0;
  end
  nMarked = nMarked + 1;
end
if ~isempty(missing)
  warning('kalmanrates:CleanMissing', ...
          '%d flagged quotes have no row in EUR_raw.mat and are not removed: %s', ...
          numel(missing), strjoin(missing, ', '));
end
save(dataPath('EUR_clean.mat'), '-struct', 'clean');
fprintf('\nSaved %s: %d quotes removed\n', dataPath('EUR_clean.mat'), nMarked);

% Curves built from an earlier cleaning no longer match
if isfile(dataPath('EUR100_clean.mat'))
  movefile(dataPath('EUR100_clean.mat'), dataPath('EUR100_clean_stale.mat'), 'f');
  fprintf('Moved the earlier EUR100_clean.mat to EUR100_clean_stale.mat (built from older flags)\n');
end

settings = struct('threshold', threshold, 'maxFlagFraction', maxFlagFraction, ...
                  'resultFile', resultFile, 'modelFile', R.modelFile, 'resultSavedAt', R.savedAt, ...
                  'firstDay', dayTimes(1), 'lastDay', dayTimes(end), 'nObs', nObs, 'nNew', nNew, ...
                  'createdAt', datestr(now, 'yyyy-mm-dd HH:MM'));
save(flagFile, 'flags', 'settings');
fprintf('Saved the flags in %s\n', flagFile);

% Keep the pass-1 result for comparePasses (kalmanOptimizeMain overwrites it)
pass1File = dataPath('kalmanOptimizeResult_pass1.mat');
if ~isfile(pass1File)
  copyfile(dataPath(resultFile), pass1File);
  fprintf('Saved the pass-1 result as %s\n', pass1File);
else
  fprintf('%s already exists and is kept (the first, raw pass).\n', pass1File);
end

useDataSet('clean');

fprintf(['\nNext (pass 2):\n' ...
         '  runTenorMultiIrCurves   %% optional, curves of the cleaned data (slow)\n' ...
         '  buildKalmanModel\n' ...
         '  kalmanOptimizeMain\n' ...
         '  comparePasses\n']);
