% buildKalmanModel.m
% Builds the Kalman model for the periods in kalmanOptimizeMain.config so
% that the out-of-sample test is out of sample also for the model itself:
% PCA (Ef, ef) and start parameters use only the in-sample days, while the
% model holds every day from insample_start to test_end (or the last day of
% the data) so that kalmanOptimizeMain can filter through the test period.
%
% Saves K and KI in data/kalmanModel_<insample_start>_<insample_end>_<test_end or end>.mat.
% Run it once per in-sample/test choice, then run kalmanOptimizeMain.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));

buildPeriods = experimentPeriods(readConfig( ...
    fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config')));
if ~buildPeriods.hasInsample
  error('kalmanrates:InvalidPeriods', ...
        'Set insample_start and insample_end in kalmanOptimizeMain.config to build a model.');
end

% runReducedIrCurvesInit uses buildPeriods to choose the days and the PCA cutoff.
runReducedIrCurvesInit;
runReducedIrCurvesCls;

KI.insampleStart = buildPeriods.insampleStart;
KI.insampleEnd = buildPeriods.insampleEnd;
KI.pcaCutoff = outOfSampleStartDate;
KI.buildEnd = buildPeriods.buildEnd;

modelPath = dataPath(buildPeriods.modelFile);
save(modelPath, 'K', 'KI', '-v7.3');
fprintf('Saved %s: %d days, %s to %s, PCA on days before %s\n', modelPath, K.nt, ...
        datestr(floor(KI.times(1)), 'yyyy-mm-dd'), datestr(floor(KI.times(end)), 'yyyy-mm-dd'), ...
        datestr(KI.pcaCutoff, 'yyyy-mm-dd'));
