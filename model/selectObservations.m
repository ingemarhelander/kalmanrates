function [indices, outputSuffix] = selectObservations(nAvailable, maxObservations, times)
% Select a prefix without changing observation frequency or order.
if ~(isnumeric(maxObservations) && isreal(maxObservations) && ...
     isscalar(maxObservations) && maxObservations > 0 && ...
     (maxObservations == Inf || ...
      (isfinite(maxObservations) && maxObservations == fix(maxObservations))))
  error('kalmanrates:InvalidObservationLimit', ...
        'maxObservations must be a positive integer or Inf.');
end
if nAvailable < 1
  error('kalmanrates:NoObservations', 'No observations are available.');
end
if nargin < 3
  times = [];
end
if ~isempty(times) && numel(times) ~= nAvailable
  error('kalmanrates:DateAlignment', 'Dates must match the observation count.');
end

nSelected = min(nAvailable, maxObservations);
indices = 1:nSelected;
outputSuffix = '';
if nSelected < nAvailable
  outputSuffix = sprintf('_sample%d', nSelected);
end
fprintf('Using %d of %d observations', nSelected, nAvailable);
if ~isempty(times)
  fprintf(': %s to %s', datestr(floor(times(1)), 'yyyy-mm-dd'), ...
          datestr(floor(times(nSelected)), 'yyyy-mm-dd'));
else
  fprintf(' (date metadata unavailable)');
end
fprintf('\n');
end
