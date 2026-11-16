function [K, times, ricAll] = loadKalmanModel()
% loadKalmanModel - Load the whole saved Kalman model with its dates and contracts.
%
%   [K, times, ricAll] = loadKalmanModel()
%
%   Unlike loadKalmanSample, no observations are dropped here. Select days
%   and contracts afterwards with selectDays, selectInstruments and
%   subsetKalmanModel.

saved = load(dataPath('kalmanDebug.mat'), 'K', 'KI');
if ~isfield(saved, 'K')
  saved = load(dataPath('kalmanDebugR.mat'), 'K', 'KI');
end
if ~isfield(saved, 'KI') || ~isfield(saved.KI, 'times') || ~isfield(saved.KI, 'ricAll')
  error('kalmanrates:MissingMetadata', ...
        'The saved model has no KI.times/KI.ricAll, so days and contracts cannot be selected.');
end
K = saved.K;
times = saved.KI.times(:);
ricAll = saved.KI.ricAll;
if numel(times) ~= K.nt
  error('kalmanrates:DateAlignment', 'KI.times must have one date per observation.');
end
end
