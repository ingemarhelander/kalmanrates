% load kalmanDebug100.mat;
load kalmanDebug.mat;
if (~exist('K', 'var')) % Only load data once (improves speed)
  load kalmanDebug.mat;
  load kalmanDebugR.mat;
end


maxIter = 100;
mu = 1E-4;
% mu = 1E-8
mu = 1E-2

thetaFInd = (1:length(K.thetaF));
SigmavInd = (1:length(K.Sigmav))+thetaFInd(end);
SigmawInd = (1:length(K.Sigmaw))+SigmavInd(end);
axInd = (1:length(K.ax))+SigmawInd(end);
Sigma0Ind = (1:length(K.Sigma0))+axInd(end);
EfVecInd = (1:length(K.Ef(:)))+Sigma0Ind(end);

load kalmanXoptAll.mat;
x = xOptAll(:,1);
K.thetaF = x(thetaFInd);
K.Sigmaw = x(SigmawInd);
K.Sigmav = x(SigmavInd);
K.ax = x(axInd);
K.Sigma0 = x(Sigma0Ind);
K.Ef(:) = x(EfVecInd);

optCase = 8;
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
end



xInitAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
xLbdAll = [-ones(size(K.thetaF))*inf ; zeros(size(K.Sigmav)) ; zeros(size(K.Sigmaw)) ; -inf*ones(size(K.ax)) ; zeros(size(K.Sigma0)) ; -inf*ones(size(K.Ef(:)))];

xLbd = xLbdAll(varInd);

fvals = [];
normOfgrads = [];
normOftravelledDistances = [];
stepSizes = [];

if(true)
    extraOutputBool = true;
    xs = [];
    grads = [];
else
    extraOutputBool = false;
end

x = xInitAll;
x = x(varInd);
xInit = x;
% x = xOpt;
fStartEnd = [];
directionalDerivativeStartEnd = [];
gNormStartEnd = [];
iter = 1;

[fval, grad] = kalmanObj(x, K, mu, varInd);

% K.thetaF = x(thetaFInd);
% K.Sigmaw = x(SigmawInd);
% K.Sigmav = x(SigmavInd);
% K.ax = x(axInd);
% K.Sigma0 = x(Sigma0Ind);
% K.Ef(:) = x(EfVecInd);
% fval2 = K.logLikelihoodEps();


while(true)
    
    fvals = [fvals;fval];
    normOfGrad = norm(grad);
    normOfgrads = [normOfgrads; normOfGrad];
    normOftravelledDistances = [normOftravelledDistances; norm(xInit - x)];
    
    if(extraOutputBool)
        grads = [grads,grad];
        xs = [xs, x];
    end

    p = grad;
    
    df = grad'*p;
    
    if (df < 0)
      fprintf('This is not an ascent direction\n');
    end
    
    alpha = 1;
    alphaLbd = ones(size(x))*inf;
    alphaUbd = ones(size(x))*inf;
    for i=1:length(p)
      if (p(i) < 0 && xLbd(i) ~= -inf)
        alphaLbd(i) = -(x(i)-xLbd(i))/p(i);
      end
    end
    
    alphaLbdMin = min(alphaLbd);
    if (alphaLbdMin < 1)
      p = p*alphaLbdMin*0.99; % Ensures that alpha = 1 is feasible
    end
    
%     alpha = min([alpha ; alphaLbdMin*0.99]);    
    
    if (true) % Steps are far to long in the beginning, sufficient with a specific increase in objective function
      objInc = 1E2;
      alphaObjInc = objInc/(grad'*p);
      alpha = min(alpha, alphaObjInc);
      alpha = max(alpha, 1E-5);
    end
    x_ = x+alpha*p;

    [fvalNew, gradNew] = kalmanObj(x_, K, mu, varInd);
    dfNew = gradNew'*p;
% K.thetaF = x_(thetaFInd);
% K.Sigmaw = x_(SigmawInd);
% K.Sigmav = x_(SigmavInd);
% K.ax = x_(axInd);
% K.Sigma0 = x_(Sigma0Ind);
% K.Ef(:) = x_(EfVecInd);
% fvalNew2 = K.logLikelihoodEps();

%     % Determine single 1-dimensional spline
%     ak = fval;
%     bk = df;
%     dxk = alpha;
%     dk = 2/dxk^3*(fval - fvalNew + (dfNew + df)/2*dxk);
%     ck = (dfNew-df)/(2*dxk) - 3*dxk/2*dk;
%     
%     tmp = ck^2-3*dk*bk;
%     if (tmp < 0 || dk == 0)
%       error('Need to investigate further');
%     end
%     alpha1 = (-ck+sqrt(tmp))/(3*dk);
%     alpha2 = (-ck-sqrt(tmp))/(3*dk);
%     if (2*ck+6*dk*alpha1 < 0)
%       alphaSpline = alpha1;
%     elseif (2*ck+6*dk*alpha2 < 0)
%       alphaSpline = alpha2;
%     else
%       error('Need to investigate further');      
%     end
%     
%     alphaMin = min([alpha1 alpha2 0]);
%     alphaMax = max([alpha1 alpha2 alpha]);
%     
%     figure(1);
%     alphaGrid = alphaMin + (alphaMax-alphaMin)*(0:0.01:1);
%     hold off;
%     plot(alphaGrid, fval + alphaGrid*df);
%     hold on;
%     plot(alphaGrid, fvalNew + (alphaGrid-alpha)*dfNew);
%     plot(alphaGrid, ak + bk*alphaGrid + ck*alphaGrid.^2 + dk*alphaGrid.^3);

    
    fLineSearch = fvalNew;
    directionalDerivativeLineSearch = gradNew'*p;
    lineSearchIter = 1;
    
fprintf('LS %2s: %12s %12s %12s\n', 'i', 'alpha', 'fvalNew-fval', 'dfNew');
% fprintf('LS %2s: %12s %12s %12s %12s %12s\n', 'i', 'alpha', 'fvalNew-fval', 'fvalNew2-fval2', 'dfNew', 'fvalNew2-fvalNew');
    while (true)
      isDescent = directionalDerivativeLineSearch(end) < 0;
%       isBetter = fvalNew > fval*(1+1E-8);
%       isBetter = fvalNew > (fval+abs(fval)*1E-8);
      isBetter = fvalNew > fval;
      dfNew = gradNew'*p;
fprintf('LS %2d: %12d %12f %12d\n', lineSearchIter, alpha, fvalNew-fval, dfNew);
% fprintf('LS %2d: %12d %12d %12f %12d %12d\n', lineSearchIter, alpha, fvalNew-fval, fvalNew2-fval2, dfNew, fvalNew2-fvalNew);
if (isnan(fvalNew))
  hej = 1
end
%       if (isDescent && isBetter && isAccurate)
%       if (isBetter && isAccurate)
      if (isBetter)
        break;
      end
      alpha = alpha/2;

      x_ = x+alpha*p;

      [fvalNew, gradNew] = kalmanObj(x_, K, mu, varInd);
% K.thetaF = x_(thetaFInd);
% K.Sigmaw = x_(SigmawInd);
% K.Sigmav = x_(SigmavInd);
% K.ax = x_(axInd);
% K.Sigma0 = x_(Sigma0Ind);
% K.Ef(:) = x_(EfVecInd);
% fvalNew2 = K.logLikelihoodEps();
      
      directionalDerivativeLineSearch = [directionalDerivativeLineSearch gradNew'*p];
      
      lineSearchIter = lineSearchIter+1;
      if (lineSearchIter > 50)
        break;
%         error('Could not find improvement');
      end
    end
    
    
   
    if (false)
%       figure(1);
%       plotObj(func, x, p*alpha, B);

      lambdaStep = ((-9:1:9)/9)'*alpha;
      lambdaStep = ((0:1:9)/9)'*alpha;
%       lambdaStep = (0:0.1:1)'*alpha;
%       lambdaStep = (0:0.1:2)'*alpha;
      fAll = zeros(size(lambdaStep));
      fApprox = zeros(size(lambdaStep));
      gNorm = zeros(size(lambdaStep));
      directionalDerivativeAll = zeros(size(lambdaStep));
      for i=1:length(lambdaStep)
        [fAll(i), gTmp] = kalmanObj(x+p*lambdaStep(i), K, mu, varInd);
       
  %       fApprox(i) = fval + gNum*lambda(i) + 0.5*hNum*lambda(i)^2;
        dx = p*lambdaStep(i);
        fApprox(i) = fval + grad'*dx + 0.5*dx'*B*dx;          
        gNorm(i) = norm(gTmp);
        directionalDerivativeAll(i) = gTmp'*p;
      end
      fStartEnd = [fStartEnd ; fAll];
      gNormStartEnd = [gNormStartEnd ; gNorm];
      directionalDerivativeStartEnd = [directionalDerivativeStartEnd ; directionalDerivativeAll];
      figure(1);
      clf
      plot(fStartEnd)
      figure(2);
      clf
      semilogy(gNormStartEnd)
      figure(3);
      clf
      plot(directionalDerivativeStartEnd)
    end
%     clf
%     yyaxis left;
%     plot(lambda, [fAll fApprox], alpha, fvalNew, 'g*');
%     yyaxis right;
%     plot(lambda, gNorm);
    
    if (sum(isnan(p))>0)
      stop
    end
    s = alpha*p;
    
    x = x + s;   
    
    fprintf('%3d %14d %14f %14f %14d %14d\n', iter, alpha, fvalNew-fval, fvalNew, norm(grad), norm(gradNew));
    
    fval = fvalNew;
    grad = gradNew;
    iter = iter+1;
% fval2 = fvalNew2;
    

%     x = x + dx*adjustedSteplength;
%     stepSizes = [stepSizes; adjustedSteplength];
    stepSizes = [stepSizes; alpha];

    if(norm(gradNew) < 1e-7)
        fprintf("The norm of the gradient: %+.5e\n", norm(gradL))
        break;
    end
    if(iter >= maxIter)
        fprintf("Maximum number of iterations\n")
        break;
    end
end

