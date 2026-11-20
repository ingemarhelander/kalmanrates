function KS = subsetKalmanModel(K, indices)
% subsetKalmanModel - Copy of a Kalman model restricted to consecutive days.
%
%   KS = subsetKalmanModel(K, indices)
%
%   K        clsKalmanIr object (or a struct made by this function).
%   indices  consecutive observation indices to keep (from selectDays).
%
%   KS is a struct with the same fields as K, so K itself is never changed
%   and several subsets can be made from one loaded model. kalmanObj and
%   the optimizers accept it in place of K. Parameters and states keep
%   their sizes.

indices = indices(:)';
if isempty(indices) || any(diff(indices) ~= 1) || indices(1) < 1 || indices(end) > K.nt
  error('kalmanrates:InvalidIndices', 'indices must be consecutive and within 1..%d.', K.nt);
end

if isobject(K)
  names = properties(K);
else
  names = fieldnames(K);
end
KS = struct();
for i = 1:numel(names)
  KS.(names{i}) = K.(names{i});
end

timeFields = {'x', 'w', 'z', 'v', 'F', 'Iz', 'Cw', 'Cx', 'Cv', 'nuv', ...
              'A', 'B', 'D', 'G', 'atAll', 'oAll', 'oTAll', 'oIndAll', ...
              'tcAll', 'datesStepRelAll', 'nz'};
for i = 1:numel(timeFields)
  KS.(timeFields{i}) = K.(timeFields{i})(indices);
end
KS.nt = numel(indices);
end
