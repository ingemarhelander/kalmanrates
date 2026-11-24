function r = filterDays(KS, paramVec)
% filterDays - Kalman filter as clsKalmanIr.logLikelihoodEps, for the struct
% that subsetKalmanModel returns, with the parameters in paramVec.
%
%   Per day: log likelihood f, number of observations nz, standardized squared
%   error d2 = eps'*inv(S)*eps/nz, and per contract the observed z and the
%   one-step-ahead prediction pred (NaN where the contract has no quote).

  thetaFInd = (1:length(KS.thetaF));
  SigmavInd = (1:length(KS.Sigmav)) + thetaFInd(end);
  SigmawInd = (1:length(KS.Sigmaw)) + SigmavInd(end);
  axInd     = (1:length(KS.ax))     + SigmawInd(end);
  Sigma0Ind = (1:length(KS.Sigma0)) + axInd(end);
  EfVecInd  = (1:length(KS.Ef(:)))  + Sigma0Ind(end);
  thetaF = paramVec(thetaFInd);
  Sigmav = paramVec(SigmavInd);
  Sigmaw = paramVec(SigmawInd);
  ax0    = paramVec(axInd);
  Sigma0 = paramVec(Sigma0Ind);
  Ef = KS.Ef; Ef(:) = paramVec(EfVecInd);
  EfInt = [zeros(1, size(Ef,2)) ; cumsum(Ef,1)]/365;

  nt = KS.nt;
  r.f = zeros(nt,1); r.nz = zeros(nt,1); r.d2 = zeros(nt,1);
  r.z = nan(nt, KS.nu); r.pred = nan(nt, KS.nu);
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
    S = H*PPred*H' + Cv;
    Kg = (PPred*H')/S;
    xUpdate = xPred + Kg*eps;
    PUpdate = PPred - Kg*H*PPred;

    w = S \ eps;
    r.nz(t) = KS.nz(t);
    r.f(t) = -0.5*(KS.nz(t)*log(2*pi) + log(det(S)) + eps'*w);
    r.d2(t) = (eps'*w) / KS.nz(t);
    [row, col] = find(KS.Iz{t});
    r.z(t, col) = KS.z{t}(row);
    r.pred(t, col) = zPred(row);
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
