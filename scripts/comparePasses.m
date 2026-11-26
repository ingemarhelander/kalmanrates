% comparePasses.m
% Checks the cleaning and compares pass 1 and pass 2 on the same quotes.
%
% Run after cleanMarketData, buildKalmanModel and kalmanOptimizeMain on the
% cleaned data (useDataSet('clean')).
%
% 1. Checks that the quotes flagged by cleanMarketData (data/noisyQuotes.mat)
%    are missing in the model that is now built.
% 2. Evaluates the parameters of pass 1 (data/kalmanOptimizeResult_pass1.mat)
%    and pass 2 (data/kalmanOptimizeResult.mat) on the same model, the
%    cleaned one, so that the log likelihoods are computed on the same quotes
%    and can be compared. The log likelihood that kalmanOptimizeMain prints
%    in pass 1 and pass 2 cannot be compared, as it is computed on
%    different quotes.
%
% Without runTenorMultiIrCurves on the cleaned data the two models differ
% only by the removed quotes and the comparison is exact. With new curves
% (EUR100_clean.mat) the PCA of the model has changed as well, and the
% pass-1 parameters are used with the new model; the comparison then also
% includes the effect of the new PCA.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));

pass1File = dataPath('kalmanOptimizeResult_pass1.mat');
pass2File = dataPath('kalmanOptimizeResult.mat');
if ~isfile(pass1File)
  error('kalmanrates:MissingResult', '%s does not exist. Run cleanMarketData first.', pass1File);
end
R1 = load(pass1File);
R2 = load(pass2File);
if strcmp(R1.savedAt, R2.savedAt) && isequal(R1.xOpt, R2.xOpt)
  error('kalmanrates:MissingResult', ...
        'kalmanOptimizeResult.mat is still the pass-1 result. Run buildKalmanModel and kalmanOptimizeMain first.');
end
if ~strcmp(R1.modelFile, R2.modelFile)
  error('kalmanrates:InvalidPeriods', 'Pass 1 used %s and pass 2 %s. Use the same config in both passes.', ...
        R1.modelFile, R2.modelFile);
end
fprintf('Active data set: %s\n', useDataSet());
fprintf('Pass 1: saved %s\nPass 2: saved %s\nModel:  %s\n', R1.savedAt, R2.savedAt, R2.modelFile);

[K, times, ricAll, KI] = loadKalmanModel(R2.modelFile);
modelDays = floor(times);

%% 1. Are the flagged quotes gone?
Q = load(dataPath('noisyQuotes.mat'), 'flags', 'settings');
nIn = 0; nPresent = 0; nNoContract = 0;
present = {};
for i = 1:height(Q.flags)
  t = find(modelDays == Q.flags.date(i), 1);
  if isempty(t), continue; end
  nIn = nIn + 1;
  c = find(strcmp(ricAll, Q.flags.ric{i}), 1);
  if isempty(c)
    nNoContract = nNoContract + 1;
    continue;
  end
  if any(K.Iz{t}(:, c))
    nPresent = nPresent + 1;
    present{end+1} = sprintf('%s %s', datestr(Q.flags.date(i), 'yyyy-mm-dd'), Q.flags.ric{i}); %#ok<SAGROW>
  end
end
fprintf('\n=== Check of the cleaning ===\n');
fprintf('%d flagged quotes on the days of the model, %d still in the model', nIn, nPresent);
if nNoContract > 0
  fprintf(' (%d belong to contracts that are no longer in the model)', nNoContract);
end
fprintf('\n');
if nPresent > 0
  warning('kalmanrates:CleanNotApplied', ...
          ['%d flagged quotes are still in the model: %s. Check that useDataSet(''clean'') is active ' ...
           'and that buildKalmanModel was run after cleanMarketData.'], ...
          nPresent, strjoin(present(1:min(10, end)), ', '));
else
  fprintf('OK: all flagged quotes are missing in the model.\n');
end

%% 2. Pass 1 and pass 2 on the same quotes
nParam = numel([K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)]);
if numel(R1.xOpt) ~= nParam || numel(R2.xOpt) ~= nParam
  error('kalmanrates:ModelAlignment', ...
        'The model has %d parameters, pass 1 %d and pass 2 %d. Have contracts disappeared?', ...
        nParam, numel(R1.xOpt), numel(R2.xOpt));
end
periods = R2.periods;
if isempty(periods.outOfSampleStart)
  fprintf('\nNo out-of-sample period in the config, only the optimization days are compared.\n');
  idx = selectDays(times, periods.optimizeStart, periods.optimizeEnd);
  KS = subsetKalmanModel(K, idx);
  f1 = kalmanObj(R1.xOpt, KS, 0, 1:nParam);
  f2 = kalmanObj(R2.xOpt, KS, 0, 1:nParam);
  nObs = sum(KS.nz);
  fprintf('%-28s %16s %16s\n', '', 'pass 1', 'pass 2');
  fprintf('%-28s %16.4f %16.4f\n', 'Optimization, log lik', f1, f2);
  fprintf('%-28s %16.6f %16.6f\n', 'Optimization, per quote', f1/nObs, f2/nObs);
else
  e1 = evaluateOutOfSample(K, times, R1.xOpt, 1:nParam, periods);
  e2 = evaluateOutOfSample(K, times, R2.xOpt, 1:nParam, periods);
  fprintf('\n=== Same (cleaned) quotes, pass-1 and pass-2 parameters ===\n');
  fprintf('%-30s %16s %16s %12s\n', '', 'pass 1', 'pass 2', 'pass 2 - 1');
  row = @(name, a, b, fmt) fprintf(['%-30s ' fmt ' ' fmt ' ' fmt '\n'], name, a, b, b - a);
  row('Out-of-sample, log lik', e1.logLik, e2.logLik, '%16.4f');
  row('Out-of-sample, per day', e1.perDay, e2.perDay, '%16.4f');
  row('Out-of-sample, per quote', e1.perObs, e2.perObs, '%16.6f');
  row('Optimization, log lik', e1.optimizeLogLik, e2.optimizeLogLik, '%16.4f');
  row('Optimization, per quote', e1.optimizePerObs, e2.optimizePerObs, '%16.6f');
  fprintf('Out-of-sample %s to %s, %d days. Larger is better.\n', ...
          datestr(e1.outOfSampleStart, 'yyyy-mm-dd'), datestr(e1.outOfSampleEnd, 'yyyy-mm-dd'), e1.nDays);
end

% Measurement noise: should be smaller after the cleaning
SigmavInd = (1:numel(K.Sigmav)) + numel(K.thetaF);
v1 = R1.xOpt(SigmavInd); v2 = R2.xOpt(SigmavInd);
ok = v1 > 0 & v2 > 0;
fprintf('\nSigmav (measurement noise), pass 2 / pass 1: median %.3f, min %.3f, max %.3f\n', ...
        median(v2(ok)./v1(ok)), min(v2(ok)./v1(ok)), max(v2(ok)./v1(ok)));
