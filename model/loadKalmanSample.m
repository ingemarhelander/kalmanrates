function [K, outputSuffix] = loadKalmanSample()
% Always load the original saved model; reduce its history in memory only.
config = modelRunConfig();
saved = load(dataPath('kalmanDebug.mat'), 'K', 'KI');
if ~isfield(saved, 'K')
  saved = load(dataPath('kalmanDebugR.mat'), 'K', 'KI');
end
times = [];
if isfield(saved, 'KI') && isfield(saved.KI, 'times')
  times = saved.KI.times;
end
[K, outputSuffix] = limitKalmanObservations(saved.K, config.maxObservations, times);
end
