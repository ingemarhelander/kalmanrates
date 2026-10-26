np = size(Ef, 2);
nc = nSteps;
nu = sum(usedInstr);
nx = np + nc + nu;

if (flgTick)
  kx = [1E5 ; 1E0 ; 1E2 ; 1E1 ; 1E-1];
  
  thetaF = ones(nx,1);
  Sigmaw = [ef*kx(1) ; ones(nc,1) ; ones(nu,1)*kx(2)];
  Sigmav = ones(nu,1)*kx(4);
  ax = zeros(nx,1);      % E[x0]
  Sigma0 = ones(nx,1);   % Cov[x0] - denoted Sigmax by Carl & Axel - should be changed to Sigma0
  Sigmax = [zeros(np,1) ; ones(nc,1)*kx(5) ; ones(nu,1)*kx(3)];
else
  kx = [1E5 ; 1E0 ; 1E2 ; 1E1];
  
  thetaF = ones(nx,1);
  Sigmaw = [ef*kx(1) ; ones(nc,1) ; ones(nu,1)*kx(2)];
  Sigmav = ones(nu,1)*kx(4);
  ax = zeros(nx,1);      % E[x0]
  Sigma0 = ones(nx,1);   % Cov[x0] - denoted Sigmax by Carl & Axel - should be changed to Sigma0
  Sigmax = [zeros(np,1) ; zeros(nc,1) ; ones(nu,1)*kx(3)];
end

if (false)
  EfVec = Ef(:);
  load kalmanDebugXopt2.mat;
  x = xOpt2;
  thetaFInd = (1:length(K.thetaF));
  SigmavInd = (1:length(K.Sigmav))+thetaFInd(end);
  SigmawInd = (1:length(K.Sigmaw))+SigmavInd(end);
  axInd = (1:length(K.ax))+SigmawInd(end);
  Sigma0Ind = (1:length(K.Sigma0))+axInd(end);
  EfVecInd = (1:length(EfVec))+Sigma0Ind(end);

  K.thetaF = x(thetaFInd);
  K.Sigmaw = x(SigmawInd);
  K.Sigmav = x(SigmavInd);
  K.ax = x(axInd);
  K.Sigma0 = x(Sigma0Ind);
end

if (false) % Reduce size
  nt = 100;
  times = times(1:nt);
  indInstrAll = indInstrAll(1:nt);
  priceAll = priceAll(1:nt);
  atAll = atAll(1:nt);
  oAll = oAll(1:nt);
  oTAll = oTAll(1:nt);
  oIndAll = oIndAll(1:nt);
  tcAll = tcAll(1:nt);
  datesStepRelAll = datesStepRelAll(1:nt);
end

if (false) % Reduce size using last data points
  nt = 1000;
  times = times(end-nt:end);
  indInstrAll = indInstrAll(end-nt:end);
  priceAll = priceAll(end-nt:end);
  atAll = atAll(end-nt:end);
  oAll = oAll(end-nt:end);
  oTAll = oTAll(end-nt:end);
  oIndAll = oIndAll(end-nt:end);
  tcAll = tcAll(end-nt:end);
  datesStepRelAll = datesStepRelAll(end-nt:end);
end

K = clsKalmanIr(times, xpInit, indInstrAll, usedInstr, priceAll, pl, atAll, oAll, oTAll, oIndAll, tcAll, nSteps, datesStepRelAll, cbEffectiveDates, thetaF, Sigmaw, Sigmav, Sigmax, ax, Sigma0, Ef);
% save('kalmanDebug', 'K', 'KI');
