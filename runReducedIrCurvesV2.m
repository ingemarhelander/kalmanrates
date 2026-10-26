% Ricatti recursion - Solves optimization problem for non-linear Rauch-Tung-Striebel smoother
% Also generates data (K, KI) that is used for other Kalman filters.

runReducedIrCurvesInit;

runReducedIrCurvesCls;


% [a,b] = mexKalmanLogLikelihood(K.z, K.A, K.B, K.D, K.G, K.Iz, K.thetaF, K.Sigmav, K.Sigmaw, K.ax, K.Sigma0, K.atAll, K.oTAll, K.oIndAll, K.tcAll, Ef, K.datesStepRelAll);
[f, obj, gradObjEf, gradObjthetaF, gradObjSigmav, gradObjSigmaw, gradObjax, gradObjSigma0, b] = mexKalmanLogLikelihood(K.z, K.A, K.B, K.D, K.G, K.Iz, K.thetaF, K.Sigmav, K.Sigmaw, K.ax, diag(K.Sigma0), K.atAll, K.oTAll, K.oIndAll, K.tcAll, K.Ef, K.datesStepRelAll);
for iter=1:5

tic
  [KR, f0, fInit] = kalmanSub(K);
toc
  
  if (iter <= 100)
    lambda = 1;
  elseif (iter < 10)
    lambda = 0.01;  
  else
    lambda = 1;  
  end
  K.x0 = K.x0 + lambda*KR.dx0;
  for t=1:K.nt
    K.x{t} = K.x{t} + lambda*KR.dx{t};
    K.w{t} = K.w{t} + lambda*KR.dw{t};
    K.v{t} = K.v{t} + lambda*KR.dv{t};

%     if (t>1)
%       diff(t) = norm(K.F{t}*K.x{t-1} + K.w{t} - K.x{t});
%     end
  end
tic
  fNew = K.evalObj();
toc
  fprintf('%.16d %.16d %d\n', f0, fNew, fNew - f0);
  if (abs(f0-fInit)< 1E-8)
    break;
  end

end

axRtsOpt = K.x{1};
save('axRtsOpt', 'axRtsOpt')

%% Display solution

xp = zeros(K.nt, K.np);
xc = zeros(K.nt, K.nc);
xu = zeros(K.nt, K.nu);
v = zeros(K.nt, K.nu);
w = zeros(K.nt, K.nx);
z = nan(K.nt, K.nx);
for t=1:K.nt
  xp(t,:) = K.x{t}(1:K.np);
  xc(t,:) = K.x{t}(K.np + (1:K.nc));
  xu(t,:) = K.x{t}(K.ns+1:end);
  v(t,indInstrAll{t}) = K.v{t};
  w(t,:) = K.w{t};
  z(t,indInstrAll{t}) = K.z{t};
end

%% Figure

figure(1);
plot(xp(:,1:3));
title('x PC');

figure(2);
plot(xc(:,1:3));
title('x step');

figure(3);
plot(xu);
title('x unsystematic');

figure(4);
plot(v);
title('v');


ind = find(~isnan(cb{:,7}));
fInc = zeros(size(times));
cbRate = zeros(size(times));
for i=2:length(ind)
  fInc(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = cb{ind(i),7}-cb{ind(i-1),7};
  cbRate(times >= datenum(cb{ind(i-1),2}) & times < datenum(cb{ind(i),2})) = cb{ind(i),7};
end
gcf = figure(5);
% plot(times, xc(:,1)*100, times, fInc);
ind = (times >= datenum(2022,01,01));
plot(times(ind), cbRate(ind) + xc(ind,1)*100, times(ind), cbRate(ind));
datetick('x', 'yyyy');
title('x step');
legend('Market implied policy rate', 'ECB policy rate')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)]);
saveas(gcf,'ecbCutsHikes.pdf')


pause(0.01);

%% Figure movie

gcf = figure(6);
for k=1:K.nt
% for k=4764:K.nt

  tradeDate = floor(times(k));
  datesStep = cbEffectiveDates(cbEffectiveDates > tradeDate);
  datesStep = datesStep(1:nSteps);
  
  Es = zeros(size(Ef,1), nSteps);
  for i=1:nSteps
    Es(datesStep(i)-tradeDate+1:end, i) = 1;
  end
  E = [Ef Es];
  
  if (flgTick)
    plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]');
    title(datestr(times(k)));
  else
    plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]', (0:size(fH,2)-1)/365, fH(k,:));
    title(datestr(floor(times(k))));
  end
  pause(0.01);
  % pause;

  if (false)
  figure(5+mod(k,2))
  plot((0:size(E,1)-1)/365, E(:,7:end)*[xc(k,:)]');
  title(datestr(floor(times(k))));
  hold off;

  figure(5+mod(k+1,2))
  plot((0:size(E,1)-1)/365, E(:,7:end)*[xc(k,:)]');
%   title(datestr(floor(times(k))));
  hold on;
  
  pause;
  end
end
