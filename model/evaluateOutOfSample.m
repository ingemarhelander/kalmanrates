function result = evaluateOutOfSample(K, times, x, varInd, periods)
% evaluateOutOfSample - Log likelihood of the out-of-sample days with fixed parameters.
%
%   result = evaluateOutOfSample(K, times, x, varInd, periods)
%
%   K          the whole model (all days) holding the start parameters.
%   times      datenum per observation in K.
%   x, varInd  the optimized parameters, as passed to kalmanObj.
%   periods    from experimentPeriods (optimizeStart/optimizeEnd,
%              outOfSampleStart/outOfSampleEnd).
%
%   The filter is run from the first optimization day to the last
%   out-of-sample day and the log likelihood of the days before the
%   out-of-sample period is subtracted. What is left is the log likelihood
%   of the out-of-sample days given everything observed before them
%   (log p(out-of-sample | earlier days)), the usual out-of-sample measure
%   for a Kalman filter. Days between the optimization and the out-of-sample
%   period only update the state, they never fit parameters.
%
%   The mu penalty is not included, so the values are pure log likelihoods.

optimizeIdx = selectDays(times, periods.optimizeStart, periods.optimizeEnd);
outOfSampleIdx = selectDays(times, periods.outOfSampleStart, periods.outOfSampleEnd);
if outOfSampleIdx(1) <= optimizeIdx(end)
  error('kalmanrates:InvalidPeriods', 'The out-of-sample period must come after the optimization period.');
end

[fAll, nAll] = logLik(K, optimizeIdx(1):outOfSampleIdx(end), x, varInd);
[fBefore, nBefore] = logLik(K, optimizeIdx(1):outOfSampleIdx(1)-1, x, varInd);
[fOptimize, nOptimize] = logLik(K, optimizeIdx, x, varInd);

result.outOfSampleStart = floor(times(outOfSampleIdx(1)));
result.outOfSampleEnd = floor(times(outOfSampleIdx(end)));
result.logLik = fAll - fBefore;
result.nDays = numel(outOfSampleIdx);
result.perDay = result.logLik / result.nDays;
result.perObs = result.logLik / (nAll - nBefore);
result.optimizeLogLik = fOptimize;
result.optimizeDays = numel(optimizeIdx);
result.optimizePerDay = fOptimize / result.optimizeDays;
result.optimizePerObs = fOptimize / nOptimize;
end

function [f, nObs] = logLik(K, indices, x, varInd)
KS = subsetKalmanModel(K, indices);
f = kalmanObj(x, KS, 0, varInd);
nObs = sum(KS.nz);
end
