run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));

[K, runOutputSuffix] = loadKalmanSample();
thetaFInd = (1:length(K.thetaF));
SigmavInd = (1:length(K.Sigmav))+thetaFInd(end);
SigmawInd = (1:length(K.Sigmaw))+SigmavInd(end);
axInd = (1:length(K.ax))+SigmawInd(end);
Sigma0Ind = (1:length(K.Sigma0))+axInd(end);
EfVecInd = (1:length(K.Ef(:)))+Sigma0Ind(end);

load(dataPath('kalmanXoptAll.mat'), 'xOptAll');
x = xOptAll(:,end);

% Generated xOptAll(:,3)
% load axRtsOpt; % Optimal solution of ax from RTS - creates much more stable solution
% x(axInd) = axRtsOpt;
% x(thetaFInd(23:end)) = 0.99;

% x = xOptAll(:,1);
% x(thetaFInd) = 1;
% x = xInitAll;

K.thetaF = x(thetaFInd);
K.Sigmaw = x(SigmawInd);
K.Sigmav = x(SigmavInd);
K.ax = x(axInd);
K.Sigma0 = x(Sigma0Ind);
K.Ef(:) = x(EfVecInd);

optCase = 6;
if (optCase == 1)
  varInd = [thetaFInd SigmavInd SigmawInd axInd Sigma0Ind EfVecInd];
elseif (optCase == 2)
  varInd = [SigmavInd SigmawInd axInd Sigma0Ind];
elseif (optCase == 3)
  varInd = [SigmavInd SigmawInd Sigma0Ind];
elseif (optCase == 4)
  varInd = [EfVecInd];
elseif (optCase == 5)
  varInd = [SigmavInd];
elseif (optCase == 6)
  varInd = [thetaFInd SigmavInd SigmawInd Sigma0Ind];
elseif (optCase == 7)
  varInd = [thetaFInd];
elseif (optCase == 8)
  varInd = [SigmawInd];
elseif (optCase == 9)
  varInd = [thetaFInd SigmavInd SigmawInd axInd Sigma0Ind];
end

mu = 0;


xInitAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
xLbdAll = [-ones(size(K.thetaF))*inf ; zeros(size(K.Sigmav)) ; zeros(size(K.Sigmaw)) ; -inf*ones(size(K.ax)) ; zeros(size(K.Sigma0)) ; -inf*ones(size(K.Ef(:)))];

x0 = x(varInd);
lb = xLbdAll(varInd);

ub = inf*ones(size(lb));
optionvec=optimset('MaxFunEvals',2000,'Display','iter','TolX',1e-12,'TolFun',1e-7,'Algorithm','interior-point','GradObj','on');

% x = [ nu    beta0   beta1  beta2]
A   = [];
b   = [];

% Solve min  f(x)
%       s.t. Ax <= b
%            lb <= x <= ub

[xOpt,f,exitflag,output,lambda,grad,hessian] = fmincon(@(x) kalmanObjMin(x, K, mu, varInd), x0,A,b,[],[],lb,ub,[],optionvec);
xOpt2 = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)]; xOpt2(varInd) = xOpt;

% xOptAll = [xOptAll xOpt]; fOptAll = [fOptAll kalmanObj(xOptAll(:,end), K, mu)]; save(dataPath(['kalmanXoptAll' runOutputSuffix '.mat']), 'xOptAll', 'fOptAll');

% xOptAll = [xOpt2]; fOptAll = [kalmanObj(xOptAll(:,end), K, mu)]; save(dataPath(['kalmanXoptAll' runOutputSuffix '.mat']), 'xOptAll', 'fOptAll');
% xOptAll = [xOptAll xOpt2]; fOptAll = [fOptAll kalmanObj(xOptAll(:,end), K, mu)]; save(dataPath(['kalmanXoptAll' runOutputSuffix '.mat']), 'xOptAll', 'fOptAll');

