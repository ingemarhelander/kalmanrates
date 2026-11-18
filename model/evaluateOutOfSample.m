function result = evaluateOutOfSample(K, times, keepInstr, x, varInd, trainOpts, testOpts)
% evaluateOutOfSample - Log likelihood of the test days with fixed parameters.
%
%   result = evaluateOutOfSample(K, times, keepInstr, x, varInd, trainOpts, testOpts)
%
%   K          the whole model (all days) holding the start parameters.
%   times      datenum per observation in K.
%   keepInstr  contracts to observe (from selectInstruments).
%   x, varInd  the optimized parameters, as passed to kalmanObj.
%   trainOpts  selectDays options of the training data.
%   testOpts   startDate/endDate of the test period. The weekday and
%              dayStep filters of trainOpts are used for the test days too.
%
%   When the test period starts after the training period starts, the
%   filter is run from the first training day up to the last test day, and
%   the log likelihood of the days before the test period is subtracted.
%   What is left is the log likelihood of the test days given everything
%   observed before them (log p(test | earlier days)), which is the usual
%   out-of-sample measure for a Kalman filter. Days between the training
%   and the test period are only used to update the state, never to fit
%   parameters. Otherwise the filter is started at the first test day.
%
%   The mu penalty is not included, so the values are pure log likelihoods.

window = selectDays(times, struct('startDate', getOpt(testOpts, 'startDate'), ...
                                  'endDate', getOpt(testOpts, 'endDate')));
testFirst = floor(times(window(1)));
testLast = floor(times(window(end)));

trainIdx = selectDays(times, trainOpts);
trainFirst = floor(times(trainIdx(1)));
trainLast = floor(times(trainIdx(end)));

result = struct();
result.testStart = testFirst;
result.testEnd = testLast;
if testFirst <= trainLast && testLast >= trainFirst
  warning('kalmanrates:OverlappingPeriods', ...
          'The test period overlaps the training period, so it is not out of sample.');
end

filters = trainOpts;
filters.endDate = testLast;
if trainFirst < testFirst
  filters.startDate = trainFirst;
  histIdx = selectDays(times, filters);
  isTest = floor(times(histIdx)) >= testFirst;
  prefixIdx = histIdx(~isTest);
  testIdx = histIdx(isTest);
  if isempty(testIdx)
    error('kalmanrates:NoObservations', 'No test days match the selection.');
  end
  [fAll, nAll] = logLik(K, histIdx, keepInstr, x, varInd);
  [fPrefix, nPrefix] = logLik(K, prefixIdx, keepInstr, x, varInd);
  result.logLik = fAll - fPrefix;
  result.nObs = nAll - nPrefix;
  result.conditional = true;
else
  filters.startDate = testFirst;
  testIdx = selectDays(times, filters);
  [result.logLik, result.nObs] = logLik(K, testIdx, keepInstr, x, varInd);
  result.conditional = false;
end

result.nDays = numel(testIdx);
result.perDay = result.logLik / result.nDays;
result.perObs = result.logLik / result.nObs;

[result.trainLogLik, nTrainObs] = logLik(K, trainIdx, keepInstr, x, varInd);
result.trainDays = numel(trainIdx);
result.trainPerDay = result.trainLogLik / result.trainDays;
result.trainPerObs = result.trainLogLik / nTrainObs;
end

function [f, nObs] = logLik(K, indices, keepInstr, x, varInd)
KS = subsetKalmanModel(K, indices, keepInstr);
f = kalmanObj(x, KS, 0, varInd);
nObs = sum(KS.nz);
end

function value = getOpt(opts, name)
value = [];
if isfield(opts, name)
  value = opts.(name);
end
end
