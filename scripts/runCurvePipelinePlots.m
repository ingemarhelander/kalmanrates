% Build, optimize, and plot the interest-rate curve.
% Set rebuildMarketCurves to false to reuse the existing EUR100.mat curves.
% Similar setup to runReducedIrCurvesV2.m, but with additional plotting of the 
%smoothed states and reconstructed curves. Including Kalman run.


%-------------- Load Optimized Parameters ---------------------------
scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
run(fullfile(repoDir, 'setupPaths.m'));

rebuildMarketCurves = true;
showCurveMovie = true;
showStepComponentMovie = true;

if rebuildMarketCurves
  run(fullfile(scriptDir, 'runTenorMultiIrCurves.m')); %This script builds the EUR100.mat curves and saves them to dataPath('EUR100.mat').
end

run(fullfile(scriptDir, 'runReducedIrCurvesInit.m')); %As in runReducedIrCurvesV2.m, this script loads the EUR100.mat curves and builds the K and KI objects.
run(fullfile(scriptDir, 'runReducedIrCurvesCls.m')); 

optimizedParametersFile = dataPath(['kalmanOptimizeMainOutput' runOutputSuffix '.mat']);
optimizerResult = load(optimizedParametersFile, 'kalmanOptimizeMainOutput');
optimizedParameters = optimizerResult.kalmanOptimizeMainOutput.parameterVector;
expectedParameterCount = numel(K.thetaF) + numel(K.Sigmav) + ...
    numel(K.Sigmaw) + numel(K.ax) + numel(K.Sigma0) + numel(K.Ef);
K.storeParamVec(optimizedParameters);
paramVec = optimizedParameters;

% ---- RTS state-smoothin from runReducedCurvesV2.m------
for iter = 1:5
  [KR, f0, fInit] = kalmanSub(K);
  if iter <= 100
    lambda = 1;
  elseif iter < 10
    lambda = 0.01;
  else
    lambda = 1;
  end
  K.x0 = K.x0 + lambda * KR.dx0;
  for t = 1:K.nt
    K.x{t} = K.x{t} + lambda * KR.dx{t};
    K.w{t} = K.w{t} + lambda * KR.dw{t};
    K.v{t} = K.v{t} + lambda * KR.dv{t};
  end
  fNew = K.evalObj();
  fprintf('%.16d %.16d %d\n', f0, fNew, fNew - f0);
  if abs(f0 - fInit) < 1E-8
    break;
  end
end

Ef = K.Ef;
axRtsOpt = K.x{1};
save(dataPath(['axRtsOpt' runOutputSuffix '.mat']), 'axRtsOpt');

% Diagnostic methods defined in model/clsKalmanIr.m.
K.plotSolution(KI, paramVec);
K.analyze(KI, paramVec);

% State extraction and component plots copied from
% runReducedIrCurvesV2.m, lines 44-88.
xp = zeros(K.nt, K.np);
xc = zeros(K.nt, K.nc);
xu = zeros(K.nt, K.nu);
v = zeros(K.nt, K.nu);
w = zeros(K.nt, K.nx);
z = nan(K.nt, K.nx);
for t = 1:K.nt
  xp(t,:) = K.x{t}(1:K.np);
  xc(t,:) = K.x{t}(K.np + (1:K.nc));
  xu(t,:) = K.x{t}(K.ns+1:end);
  v(t,indInstrAll{t}) = K.v{t};
  w(t,:) = K.w{t};
  z(t,indInstrAll{t}) = K.z{t};
end

%% Figure
figure(21);
plot(xp(:,1:3));
title('x PC');

figure(22);
plot(xc(:,1:3));
title('x step');

figure(23);
plot(xu);
title('x unsystematic');

figure(24);
plot(v);
title('v');

% ECB policy-rate plot and PDF export copied from
% runReducedIrCurvesV2.m, lines 91-111.
ind = find(~isnan(cb{:,7}));
fInc = zeros(size(times));
cbRate = zeros(size(times));
for i = 2:length(ind)
  fInc(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = ...
      cb{ind(i),7} - cb{ind(i-1),7};
  cbRate(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = ...
      cb{ind(i),7};
end
gcf = figure(25);
ind = times >= datenum(2022,01,01);
if ~any(ind)
  ind = true(size(times));
end
plot(times(ind), cbRate(ind) + xc(ind,1)*100, times(ind), cbRate(ind));
datetick('x', 'yyyy');
title('x step');
legend('Market implied policy rate', 'ECB policy rate')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
saveas(gcf, fullfile(repoDir, 'docs', ['ecbCutsHikes' runOutputSuffix '.pdf']))

pause(0.01);

% Curve movie copied from runReducedIrCurvesV2.m, lines 117-149.
if showCurveMovie % So the user can see the evolution of the curves over time.
  gcf = figure(26);
  for k = 1:K.nt
    tradeDate = floor(times(k));
    datesStep = cbEffectiveDates(cbEffectiveDates > tradeDate);
    datesStep = datesStep(1:nSteps);

    Es = zeros(size(Ef,1), nSteps);
    for i = 1:nSteps
      Es(datesStep(i)-tradeDate+1:end, i) = 1;
    end
    E = [Ef Es];

    if (flgTick)
      plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]');
      title(datestr(times(k)));
    else
      plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]', ...
          (0:size(fH,2)-1)/365, fH(k,:));
      title(datestr(floor(times(k))));
    end
    pause(0.01);

    % Optional step-component movie from the if (false) block in
    % runReducedIrCurvesV2.m; use separate figures 7-8.
    if showStepComponentMovie
      figure(7 + mod(k,2));
      plot((0:size(E,1)-1)/365, E(:,K.np+1:end)*[xc(k,:)]');
      title(datestr(floor(times(k))));
      hold off;

      figure(7 + mod(k+1,2));
      plot((0:size(E,1)-1)/365, E(:,K.np+1:end)*[xc(k,:)]');
      hold on;
    end
  end
end
