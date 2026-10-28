% Integration check using the real data and MEX files, without saving outputs
% or running the optimization loops. Uses modelRunConfig like production runs.
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
runReducedIrCurvesInit;
runReducedIrCurvesCls;
assert(K.nt == numel(times) && size(xpInit,1) == K.nt);
assert(K.np == 6 && K.nc == 16);
assert(numel(KI.times) == K.nt);
for t=1:K.nt
  assert(all(indInstrAll{t} >= 1 & indInstrAll{t} <= K.nu));
  assert(size(K.Iz{t},1) == numel(K.z{t}));
end
[KR, f0, fInit] = kalmanSub(K);
assert(all(isfinite([f0; fInit; KR.dx0])));
K.x0 = K.x0 + KR.dx0;
for t=1:K.nt
  assert(all(isfinite([KR.dx{t}; KR.dw{t}; KR.dv{t}])));
  K.x{t} = K.x{t} + KR.dx{t};
  K.w{t} = K.w{t} + KR.dw{t};
  K.v{t} = K.v{t} + KR.dv{t};
end
assert(isfinite(K.evalObj()));

% All three optimizers use this loader. Check their initial parameter choices:
% Main/SD start at the first saved optimum; fmincon starts at the last.
K = loadKalmanSample();
savedOptima = load(dataPath('kalmanXoptAll.mat'), 'xOptAll');
ntheta = numel(K.thetaF);
nv = numel(K.Sigmav);
nw = numel(K.Sigmaw);
na = numel(K.ax);
n0 = numel(K.Sigma0);
sigmawIndices = ntheta+nv+(1:nw);
for column = [1 size(savedOptima.xOptAll,2)]
  allParameters = savedOptima.xOptAll(:,column);
  assert(numel(allParameters) == ntheta+nv+nw+na+n0+numel(K.Ef));
  K.thetaF = allParameters(1:ntheta);
  K.Sigmav = allParameters(ntheta+(1:nv));
  K.Sigmaw = allParameters(sigmawIndices);
  K.ax = allParameters(ntheta+nv+nw+(1:na));
  K.Sigma0 = allParameters(ntheta+nv+nw+na+(1:n0));
  K.Ef(:) = allParameters(ntheta+nv+nw+na+n0+1:end);
  [f, grad] = kalmanObj(K.Sigmaw, K, 1e-2, sigmawIndices);
  assert(isscalar(f) && isfinite(f) && all(isfinite(grad)));
  assert(numel(grad) == nw);
end
% Check the minimization wrapper and parameter blocks used by fmincon.
varInd = [1:ntheta+nv+nw, ntheta+nv+nw+na+(1:n0)];
[fMin, gradMin] = kalmanObjMin(allParameters(varInd), K, 0, varInd);
assert(isscalar(fMin) && isfinite(fMin) && all(isfinite(gradMin)));
assert(numel(gradMin) == numel(varInd));
fprintf('Data alignment, smoother, and optimizer objective/gradient checks passed.\n');
