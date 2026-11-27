function [KR, KIR, nRemoved] = removeQuotes(K, flag, KI)
% removeQuotes - Copy of a Kalman model without the flagged quotes.
%
%   KR = removeQuotes(K, flag)
%   [KR, KIR, nRemoved] = removeQuotes(K, flag, KI)
%
%   K     clsKalmanIr object or a struct with the same fields (from
%         subsetKalmanModel, or an earlier removeQuotes).
%   flag  K.nt x K.nu logical, true for the quotes to remove (e.g. from
%         flagQuotes). A flag where the contract has no quote is ignored.
%   KI    optional, the metadata from buildKalmanModel; its per-day fields
%         are reduced the same way.
%
%   KR is a struct with the same fields as K, so K itself is never changed.
%   On each day the rows of the removed quotes are dropped from every
%   per-quote field (z, Iz, G, atAll, Cv, v, nuv) and their rows in the
%   pricing data (oAll, oTAll, tcAll with the index oIndAll). The result is
%   the model that buildKalmanModel would build on data where the quotes are
%   missing, without building it again: kalmanObj, the optimizers,
%   subsetKalmanModel and innovationsLOO accept it in place of K. The number
%   of contracts and the parameter vector are unchanged.

  if ~isequal(size(flag), [K.nt K.nu])
    error('kalmanrates:InvalidFlags', 'flag must be %d x %d (days x contracts).', K.nt, K.nu);
  end
  flag = logical(flag);

  if isobject(K)
    names = properties(K);
  else
    names = fieldnames(K);
  end
  KR = struct();
  for i = 1:numel(names)
    KR.(names{i}) = K.(names{i});
  end
  if nargin >= 3
    KIR = KI;
  else
    KIR = [];
  end

  nRemoved = 0;
  for t = find(any(flag, 2))'
    [row, col] = find(KR.Iz{t});
    nz = size(KR.Iz{t}, 1);
    contract = zeros(nz, 1);
    contract(row) = col;
    keep = ~flag(t, contract)';
    if all(keep)
      continue;
    end
    if ~any(keep)
      error('kalmanrates:InvalidFlags', 'Day %d would have no quotes left.', t);
    end
    nRemoved = nRemoved + sum(~keep);

    KR = reduceDay(KR, t, keep, {'z', 'v', 'nuv', 'Cv', 'atAll'}, {'Iz', 'G'});
    if isfield(KR, 'nz') && numel(KR.nz) >= t
      KR.nz(t) = sum(keep);
    end
    if ~isempty(KIR)
      KIR = reduceDay(KIR, t, keep, {'indInstrAll', 'priceAll', 'maturityAll', 'atAll'}, {});
    end
  end
end

function S = reduceDay(S, t, keep, perQuote, matrixRows)
% Drops the quotes ~keep of day t from the fields of S that exist.
  nz = numel(keep);
  for f = perQuote
    if isfield(S, f{1}) && iscell(S.(f{1})) && numel(S.(f{1}){t}) == nz
      S.(f{1}){t} = S.(f{1}){t}(keep);
    end
  end
  for f = matrixRows
    if isfield(S, f{1})
      S.(f{1}){t} = S.(f{1}){t}(keep, :);
    end
  end
  % Pricing rows: quote j has rows oInd(j):oInd(j+1)-1 in oAll, oTAll, tcAll
  if isfield(S, 'oIndAll') && iscell(S.oIndAll) && numel(S.oIndAll{t}) == nz + 1
    oInd = S.oIndAll{t};
    len = diff(oInd(:));
    keepRows = repelem(keep(:), len);
    for f = {'oAll', 'oTAll', 'tcAll'}
      if isfield(S, f{1}) && iscell(S.(f{1})) && size(S.(f{1}){t}, 1) == numel(keepRows)
        S.(f{1}){t} = S.(f{1}){t}(keepRows, :);
      end
    end
    newInd = [1 ; 1 + cumsum(len(keep(:)))];
    if size(oInd, 1) > 1
      S.oIndAll{t} = cast(newInd, 'like', oInd);
    else
      S.oIndAll{t} = cast(newInd', 'like', oInd);
    end
  end
end
