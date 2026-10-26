function [f, grad] = kalmanObj(x, K, mu, varInd)

if (nargin == 3)
  varInd = 1:length(x);
end

thetaFInd = (1:length(K.thetaF));
SigmavInd = (1:length(K.Sigmav))+thetaFInd(end);
SigmawInd = (1:length(K.Sigmaw))+SigmavInd(end);
axInd = (1:length(K.ax))+SigmawInd(end);
Sigma0Ind = (1:length(K.Sigma0))+axInd(end);
EfVecInd = (1:length(K.Ef(:)))+Sigma0Ind(end);

xAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
xAll(varInd) = x;

K.thetaF = xAll(thetaFInd);
K.Sigmaw = xAll(SigmawInd);
K.Sigmav = xAll(SigmavInd);
K.ax = xAll(axInd);
K.Sigma0 = xAll(Sigma0Ind); % Currently only a diagonal matrix is being used
K.Ef(:) = xAll(EfVecInd);

if (nargout == 1)
  [f] = mexKalmanLogLikelihood(K.z, K.A, K.B, K.D, K.G, K.Iz, K.thetaF, K.Sigmav, K.Sigmaw, K.ax, diag(K.Sigma0), K.atAll, K.oTAll, K.oIndAll, K.tcAll, K.Ef, K.datesStepRelAll);

  if (mu~=0)
    f = f + mu*(sum(log(K.Sigmaw)) + sum(log(K.Sigmav)) + sum(log(K.Sigma0)));
  end
else
  [f, obj, gradObjEf, gradObjthetaF, gradObjSigmav, gradObjSigmaw, gradObjax, gradObjSigma0, b] = mexKalmanLogLikelihood(K.z, K.A, K.B, K.D, K.G, K.Iz, K.thetaF, K.Sigmav, K.Sigmaw, K.ax, diag(K.Sigma0), K.atAll, K.oTAll, K.oIndAll, K.tcAll, K.Ef, K.datesStepRelAll);

  if (mu~=0)
    f = f + mu*(sum(log(K.Sigmaw)) + sum(log(K.Sigmav)) + sum(log(K.Sigma0)));
    gradObjSigmav = gradObjSigmav + mu./K.Sigmav;
    gradObjSigmaw = gradObjSigmaw + mu./K.Sigmaw;
    gradObjSigma0 = gradObjSigma0 + mu./K.Sigma0;
  end

  grad = [gradObjthetaF ; gradObjSigmav ; gradObjSigmaw ; gradObjax ; diag(gradObjSigma0) ; gradObjEf(:)];
  grad = grad(varInd);
end

