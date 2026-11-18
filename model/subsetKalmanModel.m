function KS = subsetKalmanModel(K, indices, keepInstr)
% subsetKalmanModel - Copy of a Kalman model restricted to some days and contracts.
%
%   KS = subsetKalmanModel(K, indices, keepInstr)
%
%   K          clsKalmanIr object (or a struct made by this function).
%   indices    increasing observation indices to keep (from selectDays).
%   keepInstr  logical 1 x K.nu, the instruments to observe (from
%              selectInstruments). Default: all.
%
%   KS is a struct with the same fields as K, so K itself is never changed
%   and several subsets can be made from one loaded model. kalmanObj and
%   the optimizers accept it in place of K.
%
%   The state vector and the parameter vector keep their sizes. Dropping a
%   contract only removes its observations; its state is still propagated.
%
%   When days are skipped, a central bank step (the shift in K.B) that falls
%   on a skipped day is moved to the next kept day. At most one step may
%   fall in each gap, since the transition only shifts one step at a time.

if nargin < 3 || isempty(keepInstr)
  keepInstr = true(1, K.nu);
end
keepInstr = logical(keepInstr(:)');
if numel(keepInstr) ~= K.nu
  error('kalmanrates:ModelAlignment', 'keepInstr must have one value per instrument (%d).', K.nu);
end
if ~any(keepInstr)
  error('kalmanrates:NoInstruments', 'At least one contract must be selected.');
end
indices = indices(:)';
if isempty(indices) || any(diff(indices) <= 0) || indices(1) < 1 || indices(end) > K.nt
  error('kalmanrates:InvalidIndices', 'indices must be increasing and within 1..%d.', K.nt);
end

timeFields = {'x', 'w', 'z', 'v', 'F', 'Iz', 'Cw', 'Cx', 'Cv', 'nuv', ...
              'A', 'B', 'D', 'G', 'atAll', 'oAll', 'oTAll', 'oIndAll', ...
              'tcAll', 'datesStepRelAll'};
for i = 1:numel(timeFields)
  values = K.(timeFields{i});
  if ~iscell(values) || numel(values) ~= K.nt
    error('kalmanrates:ModelAlignment', ...
          'K.%s must contain one cell per observation.', timeFields{i});
  end
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

% Central bank steps: which days shift the step states.
identity = speye(K.nx);
isShift = false(1, K.nt);
for t = 1:K.nt
  isShift(t) = nnz(K.B{t} - identity) > 0;
end
newB = K.B(indices);
for j = 2:numel(indices)
  gap = (indices(j-1)+1):indices(j);
  nShift = sum(isShift(gap));
  if nShift > 1
    error('kalmanrates:TooManySteps', ...
          ['%d central bank steps fall between observations %d and %d. Keep ' ...
           'more days (smaller day_step) so that at most one falls between kept days.'], ...
          nShift, indices(j-1), indices(j));
  elseif nShift == 1
    newB{j} = K.B{gap(isShift(gap))};
  else
    newB{j} = identity;
  end
end

for i = 1:numel(timeFields)
  KS.(timeFields{i}) = K.(timeFields{i})(indices);
end
KS.B = newB;
KS.nz = K.nz(indices);
KS.nt = numel(indices);
for j = 1:KS.nt
  KS.F{j} = sparse(KS.A{j} * diag(KS.thetaF) * KS.B{j});
end

if all(keepInstr)
  return;
end
for j = 1:KS.nt
  KS = dropInstruments(KS, j, keepInstr);
end
end

function KS = dropInstruments(KS, t, keepInstr)
% Remove the observation rows of unselected instruments on day t.
nz = KS.nz(t);
[row, col] = find(KS.Iz{t});
instrOfRow = zeros(nz, 1);
instrOfRow(row) = col;
keepRow = keepInstr(instrOfRow)';
if ~any(keepRow)
  error('kalmanrates:NoObservations', ...
        'Observation %d has none of the selected contracts.', t);
end

% Each instrument owns a block of rows in o, oT and tc (oIndAll gives the starts).
oInd = KS.oIndAll{t}(:);
keepO = false(oInd(end)-1, 1);
blockLen = diff(oInd);
for r = find(keepRow)'
  keepO(oInd(r):oInd(r+1)-1) = true;
end

perObs = {'z', 'v', 'Cv', 'nuv', 'atAll'};
for i = 1:numel(perObs)
  value = KS.(perObs{i}){t};
  if numel(value) == nz
    KS.(perObs{i}){t} = value(keepRow);
  end
end
perRow = {'oTAll', 'tcAll', 'oAll'};
for i = 1:numel(perRow)
  value = KS.(perRow{i}){t};
  if size(value, 1) == numel(keepO)
    KS.(perRow{i}){t} = value(keepO, :);
  end
end
KS.Iz{t} = KS.Iz{t}(keepRow, :);
KS.G{t} = KS.G{t}(keepRow, :);
KS.oIndAll{t} = [1; 1 + cumsum(blockLen(keepRow))];
KS.nz(t) = sum(keepRow);
end
