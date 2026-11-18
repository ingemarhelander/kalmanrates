function result = evaluateOutOfSample(K, times, x, varInd, periods)
% evaluateOutOfSample - Log likelihood of the test days with fixed parameters.
%
%   result = evaluateOutOfSample(K, times, x, varInd, periods)
%
%   K          the whole model (all days) holding the start parameters.
%   times      datenum per observation in K.
%   x, varInd  the optimized parameters, as passed to kalmanObj.
%   periods    from experimentPeriods (trainStart/trainEnd, testStart/testEnd).
%
%   The filter is run from the first training day to the last test day and
%   the log likelihood of the days before the test period is subtracted.
%   What is left is the log likelihood of the test days given everything
%   observed before them (log p(test | earlier days)), the usual
%   out-of-sample measure for a Kalman filter. Days between the training
%   and the test period only update the state, they never fit parameters.
%
%   The mu penalty is not included, so the values are pure log likelihoods.

trainIdx = selectDays(times, periods.trainStart, periods.trainEnd);
testIdx = selectDays(times, periods.testStart, periods.testEnd);
if testIdx(1) <= trainIdx(end)
  error('kalmanrates:InvalidPeriods', 'The test period must come after the training period.');
end

[fAll, nAll] = logLik(K, trainIdx(1):testIdx(end), x, varInd);
[fBefore, nBefore] = logLik(K, trainIdx(1):testIdx(1)-1, x, varInd);
[fTrain, nTrain] = logLik(K, trainIdx, x, varInd);

result.testStart = floor(times(testIdx(1)));
result.testEnd = floor(times(testIdx(end)));
result.logLik = fAll - fBefore;
result.nDays = numel(testIdx);
result.perDay = result.logLik / result.nDays;
result.perObs = result.logLik / (nAll - nBefore);
result.trainLogLik = fTrain;
result.trainDays = numel(trainIdx);
result.trainPerDay = fTrain / result.trainDays;
result.trainPerObs = fTrain / nTrain;
end

function [f, nObs] = logLik(K, indices, x, varInd)
KS = subsetKalmanModel(K, indices);
f = kalmanObj(x, KS, 0, varInd);
nObs = sum(KS.nz);
end
