function r = innovationsLOO(KS, paramVec, opts)
% innovationsLOO - Kalman filter that flags quotes that deviate on their own.
%
%   r = innovationsLOO(KS, paramVec)
%   r = innovationsLOO(KS, paramVec, opts)
%
%   The same filter as filterDays (and clsKalmanIr.logLikelihoodEps), for a
%   model K or a struct from subsetKalmanModel with the whole parameter
%   vector paramVec (as storeParamVec takes it). The filter is not changed.
%   For each day the leave-one-out innovation of every quote is computed,
%   the error of quote i given the past AND all other quotes the same day:
%
%       w = S\eps,   rLOO_i = w_i / sqrt((S^-1)_ii)
%
%   which is N(0,1) when the model holds. When many contracts move together
%   (a rate change, a crisis day) the other quotes explain the move and
%   rLOO stays small. A quote that jumps on its own gets a large rLOO.
%
%   The worst quote with |rLOO| > threshold is removed, rLOO is computed
%   again for the remaining quotes (two bad quotes the same day can otherwise
%   mask each other), and so on. The update step of the day then uses only
%   the remaining quotes, so a removed quote does not affect later days.
%
%   opts (all optional)
%     threshold        flag when |rLOO| > threshold. Default 4.
%     maxFlagFraction  at most this share of a day's quotes is removed.
%                      Default 0.2. Days that reach it are marked in
%                      r.capped: there the model does not explain the day
%                      and the day should be looked at by hand.
%     protected        1 x nu logical, contracts that are never removed
%                      (e.g. the overnight fixing). Default none.
%
%   r has, per day t and contract i (NaN where the contract has no quote):
%     z, pred      observed quote and one-step-ahead prediction
%     rStd         ordinary standardized innovation eps_i/sqrt(S_ii)
%     rLOO         leave-one-out innovation with all quotes of the day
%     flag         logical, the quote was removed
%     rRemoved     rLOO when the quote was removed (given the quotes left
%                  at that point; can differ from rLOO when another bad
%                  quote the same day masked it)
%   and per day:
%     nz           number of quotes, nFlag removed, capped (see above)
%     f            log likelihood of the remaining quotes
%     d2           eps'*inv(S)*eps/nz of all quotes (as filterDays)

  if nargin < 3, opts = struct(); end
  threshold = getOpt(opts, 'threshold', 4);
  maxFlagFraction = getOpt(opts, 'maxFlagFraction', 0.2);
  protected = getOpt(opts, 'protected', false(1, KS.nu));
  if numel(protected) ~= KS.nu
    error('kalmanrates:InvalidOption', 'protected must have one element per contract (%d).', KS.nu);
  end
  protected = logical(protected(:)');

  thetaFInd = (1:length(KS.thetaF));
  SigmavInd = (1:length(KS.Sigmav)) + thetaFInd(end);
  SigmawInd = (1:length(KS.Sigmaw)) + SigmavInd(end);
  axInd     = (1:length(KS.ax))     + SigmawInd(end);
  Sigma0Ind = (1:length(KS.Sigma0)) + axInd(end);
  EfVecInd  = (1:length(KS.Ef(:)))  + Sigma0Ind(end);
  if numel(paramVec) ~= EfVecInd(end)
    error('kalmanrates:ModelAlignment', 'paramVec has %d elements but the model needs %d.', ...
          numel(paramVec), EfVecInd(end));
  end
  paramVec = paramVec(:);
  thetaF = paramVec(thetaFInd);
  Sigmav = paramVec(SigmavInd);
  Sigmaw = paramVec(SigmawInd);
  ax0    = paramVec(axInd);
  Sigma0 = paramVec(Sigma0Ind);
  Ef = KS.Ef; Ef(:) = paramVec(EfVecInd);
  EfInt = [zeros(1, size(Ef,2)) ; cumsum(Ef,1)]/365;

  nt = KS.nt;
  nu = KS.nu;
  r.z = nan(nt, nu); r.pred = nan(nt, nu);
  r.rStd = nan(nt, nu); r.rLOO = nan(nt, nu);
  r.flag = false(nt, nu); r.rRemoved = nan(nt, nu);
  r.nz = zeros(nt,1); r.nFlag = zeros(nt,1); r.capped = false(nt,1);
  r.f = zeros(nt,1); r.d2 = zeros(nt,1);
  r.threshold = threshold;
  r.maxFlagFraction = maxFlagFraction;

  for t = 1:nt
    F  = sparse(KS.A{t}*diag(thetaF)*KS.B{t});
    Cv = KS.G{t}*diag(Sigmav)*KS.G{t}';
    Cw = KS.D{t}*diag(Sigmaw)*KS.D{t}';
    if t == 1
      xPred = F*ax0;
      PPred = F*diag(Sigma0)*F' + Cw;
    else
      xPred = F*xUpdate;
      PPred = F*PUpdate*F' + Cw;
    end
    o = oCompute(KS, t, EfInt);
    [g, G] = taylorInstrPC(KS.pl, KS.atAll{t}, o, KS.oIndAll{t}, KS.tcAll{t}, xPred(1:KS.ns));
    u = g - G*xPred(1:KS.ns);
    H = [G KS.Iz{t}];
    zPred = H*xPred + u;
    eps = KS.z{t} - zPred;
    S = full(H*PPred*H' + Cv);
    S = (S + S')/2;

    nz = numel(eps);
    [row, col] = find(KS.Iz{t});
    contract = zeros(nz, 1);
    contract(row) = col;

    r.nz(t) = nz;
    r.z(t, contract) = KS.z{t};
    r.pred(t, contract) = zPred;
    r.rStd(t, contract) = eps ./ sqrt(diag(S));
    r.d2(t) = (eps' * (S \ eps)) / nz;

    % Remove the worst quote while it deviates more than the threshold
    keep = true(nz, 1);
    maxRemove = floor(maxFlagFraction * nz);
    first = true;
    while true
      idx = find(keep);
      Sinv = inv(S(idx, idx));
      Sinv = (Sinv + Sinv')/2;
      rl = (Sinv * eps(idx)) ./ sqrt(diag(Sinv));
      if first
        r.rLOO(t, contract) = rl;
        first = false;
      end
      cand = abs(rl);
      cand(protected(contract(idx))) = -Inf;
      [m, j] = max(cand);
      if isempty(m) || m <= threshold || numel(idx) <= 1
        break;
      end
      if nz - numel(idx) >= maxRemove
        r.capped(t) = true;
        break;
      end
      keep(idx(j)) = false;
      r.flag(t, contract(idx(j))) = true;
      r.rRemoved(t, contract(idx(j))) = rl(j);
    end
    r.nFlag(t) = sum(~keep);

    % Update with the remaining quotes only
    Hk = H(keep, :);
    Sk = S(keep, keep);
    ek = eps(keep);
    Kg = (PPred*Hk')/Sk;
    xUpdate = xPred + Kg*ek;
    PUpdate = PPred - Kg*Hk*PPred;
    r.f(t) = -0.5*(numel(ek)*log(2*pi) + log(det(Sk)) + ek'*(Sk\ek));
  end
end

function v = getOpt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end

function o = oCompute(KS, t, EfInt)
% Same as the private clsKalmanIr.oCompute.
  datesStepRel = KS.datesStepRelAll{t};
  o = [-EfInt(KS.oTAll{t},:) zeros(size(KS.oTAll{t},1), KS.nc)];
  for i = 1:KS.nc
    o(:,KS.np+i) = -max(KS.oTAll{t}-(datesStepRel(i)+1), 0)/365;
  end
end
