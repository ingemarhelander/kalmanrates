function [flag, r] = flagQuotes(K, paramVec, days, opts)
% flagQuotes - Flag quotes that deviate on their own, also with uncalibrated parameters.
%
%   flag = flagQuotes(K, paramVec, days)
%   [flag, r] = flagQuotes(K, paramVec, days, opts)
%
%   K         the whole model (clsKalmanIr or a struct from subsetKalmanModel
%             or removeQuotes).
%   paramVec  the whole parameter vector (as storeParamVec takes it), e.g.
%             the start parameters of buildKalmanModel.
%   days      consecutive day indices in K to filter (e.g. from the first
%             optimization day to the last day of the model).
%
%   The filter is run over the days with innovationsLOO. A quote is flagged
%   when its leave-one-out innovation (the error given the past and all
%   other quotes the same day) deviates more than opts.threshold times the
%   typical size for its contract.
%
%   With scale 'robust' (default) the typical size is measured from the data:
%   the filter is first run without removing anything, and the spread of
%   rLOO per contract is 1.4826*MAD (insensitive to the outliers themselves).
%   rLOO is then divided by it, so the flags do not depend on whether
%   Sigmav and Sigmaw have the right size, only on the model predicting
%   roughly right. With scale 'absolute' rLOO is taken as N(0,1), which
%   needs calibrated (optimized) parameters, as in cleanMarketData.
%
%   opts (all optional)
%     scale            'robust' (default) or 'absolute'
%     threshold        default 4
%     maxFlagFraction  default 0.2 (see innovationsLOO)
%     protected        1 x nu logical, contracts never removed. Default none.
%     burnIn           days at the start without flags. Default 20.
%     scaleDays        logical over days (numel(days)), the days used to
%                      measure the spread, e.g. only the optimization days
%                      so that later days do not affect it. Default all.
%     minQuotes        a contract with fewer quotes than this (after burnIn,
%                      on scaleDays) uses the spread of all contracts.
%                      Default 30.
%
%   flag  K.nt x K.nu logical, true for the flagged quotes (false outside days)
%   r     the result of innovationsLOO for the days (r.rScale is the spread
%         used, r.rLOO is not divided by it)

  if nargin < 4, opts = struct(); end
  scale = lower(getOpt(opts, 'scale', 'robust'));
  burnIn = getOpt(opts, 'burnIn', 20);
  minQuotes = getOpt(opts, 'minQuotes', 30);
  days = days(:)';
  nd = numel(days);
  scaleDays = logical(getOpt(opts, 'scaleDays', true(1, nd)));
  if numel(scaleDays) ~= nd
    error('kalmanrates:InvalidOption', 'scaleDays must have one element per day (%d).', nd);
  end

  KS = subsetKalmanModel(K, days);
  looOpts = struct('threshold', getOpt(opts, 'threshold', 4), ...
                   'maxFlagFraction', getOpt(opts, 'maxFlagFraction', 0.2), ...
                   'protected', getOpt(opts, 'protected', false(1, KS.nu)), ...
                   'burnIn', burnIn);

  switch scale
    case 'absolute'
      looOpts.rScale = ones(1, KS.nu);
    case 'robust'
      r0 = innovationsLOO(KS, paramVec, struct('threshold', Inf));
      rows = scaleDays(:) & ((1:nd)' > burnIn);
      R = r0.rLOO(rows, :);
      allValues = R(~isnan(R));
      if numel(allValues) < minQuotes
        error('kalmanrates:NoObservations', ...
              'Only %d quotes after burnIn to measure the spread of rLOO.', numel(allValues));
      end
      common = robustStd(allValues);
      looOpts.rScale = repmat(common, 1, KS.nu);
      for i = 1:KS.nu
        v = R(~isnan(R(:, i)), i);
        if numel(v) >= minQuotes
          looOpts.rScale(i) = robustStd(v);
        end
      end
      % A contract whose quotes are all identical to the prediction would
      % get spread 0; fall back to the common spread.
      bad = ~(looOpts.rScale > 0) | ~isfinite(looOpts.rScale);
      looOpts.rScale(bad) = common;
    otherwise
      error('kalmanrates:InvalidOption', 'scale must be ''robust'' or ''absolute''.');
  end

  r = innovationsLOO(KS, paramVec, looOpts);
  flag = false(K.nt, K.nu);
  flag(days, :) = r.flag;
end

function s = robustStd(v)
  s = 1.4826 * median(abs(v - median(v)));
end

function v = getOpt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
