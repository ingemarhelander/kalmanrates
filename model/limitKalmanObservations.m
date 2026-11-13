function [K, outputSuffix] = limitKalmanObservations(K, maxObservations, times)
% Trim only the time axis. Keep saved parameters, state dimensions and prior.
if nargin < 3
  times = [];
end
[indices, outputSuffix] = selectObservations(K.nt, maxObservations, times);
if numel(indices) == K.nt
  return;
end

timeFields = {'x', 'w', 'z', 'v', 'F', 'Iz', 'Cw', 'Cx', 'Cv', 'nuv', ...
              'A', 'B', 'D', 'G', 'atAll', 'oAll', 'oTAll', 'oIndAll', ...
              'tcAll', 'datesStepRelAll'};
% Validate before mutating the handle object, so an alignment error is atomic.
for i = 1:numel(timeFields)
  values = K.(timeFields{i});
  if ~iscell(values) || numel(values) ~= K.nt
    error('kalmanrates:ModelAlignment', ...
          'K.%s must contain one cell per observation.', timeFields{i});
  end
end
if numel(K.nz) ~= K.nt
  error('kalmanrates:ModelAlignment', 'K.nz must match K.nt.');
end
for i = 1:numel(timeFields)
  values = K.(timeFields{i});
  K.(timeFields{i}) = values(indices);
end
K.nz = K.nz(indices);
K.nt = numel(indices);
end
