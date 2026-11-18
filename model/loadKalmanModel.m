function [K, times, ricAll, KI] = loadKalmanModel(fileName)
% loadKalmanModel - Load a whole saved Kalman model with its dates and contracts.
%
%   [K, times, ricAll, KI] = loadKalmanModel()
%   [K, times, ricAll, KI] = loadKalmanModel(fileName)
%
%   fileName  file in data/, e.g. one made by buildKalmanModel. Without it
%             (or with '') kalmanDebug.mat is loaded, or kalmanDebugR.mat
%             when kalmanDebug.mat has no K.
%
%   Unlike loadKalmanSample, no observations are dropped here. Select days
%   afterwards with selectDays and subsetKalmanModel. KI is the saved
%   metadata (times, ricAll and, for
%   models from buildKalmanModel, insampleStart, insampleEnd, pcaCutoff
%   and buildEnd).

if nargin < 1 || isempty(fileName)
  saved = load(dataPath('kalmanDebug.mat'), 'K', 'KI');
  if ~isfield(saved, 'K')
    saved = load(dataPath('kalmanDebugR.mat'), 'K', 'KI');
  end
else
  path = dataPath(fileName);
  if ~isfile(path)
    error('kalmanrates:MissingModel', ...
          'The model %s does not exist. Run buildKalmanModel with the same config first.', path);
  end
  saved = load(path, 'K', 'KI');
  if ~isfield(saved, 'K')
    error('kalmanrates:MissingModel', 'The file %s has no K.', path);
  end
end
if ~isfield(saved, 'KI') || ~isfield(saved.KI, 'times') || ~isfield(saved.KI, 'ricAll')
  error('kalmanrates:MissingMetadata', ...
        'The saved model has no KI.times/KI.ricAll, so days and contracts cannot be selected.');
end
K = saved.K;
KI = saved.KI;
times = KI.times(:);
ricAll = KI.ricAll;
if numel(times) ~= K.nt
  error('kalmanrates:DateAlignment', 'KI.times must have one date per observation.');
end
end
