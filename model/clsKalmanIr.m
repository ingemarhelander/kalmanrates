% Class for Kalman filter of interest rates

classdef clsKalmanIr < handle
   properties
      nt               = 0   %
      np               = 0   % Number of principal compontents
      nc               = 0   % Number of central bank cuts/hikes
      ns               = 0   % Number of systematic risk factors (np+nc)
      nu               = 0   % Number of unsystematic risk factors
      nx               = 0   %
      nz               = 0   %

      x0               = 0   %
      x                = 0   %
      w                = 0   %
      z                = 0   %
      v                = 0   %
      F                = 0   %
      Iz               = 0   %
      Cw               = 0   %
      Cx               = 0   %
      Cv               = 0   %
      nuv              = 0   %

      thetaF           = 0   %
      Sigmaw           = 0   %
      Sigmav           = 0   %
      ax               = 0   % E[x0]
      Sigma0           = 0   % Cov[x0]
      Sigmax           = 0   % Cov[xt]

      % Matrices for keeping track of parameters
      A                = 0  % ThetaF - Rows
      B                = 0  % ThetaF - Columns
      D                = 0  % Sigmaw
      G                = 0  % Sigmav
      
      Ef               = 0   % Eigenvectors for prinicpal components
      EfInt            = 0   % Integral of Ef - used to speed up computations

      vStudentsT       = 0   %
      pl               = 0   %
      atAll            = 0   %
      oAll             = 0   %
      oTAll            = 0   %
      oIndAll          = 0   %
      tcAll            = 0   %
      datesStepRelAll  = 0
   end
   
   methods
      
        function K = clsKalmanIr(times, xpInit, indInstrAll, usedInstr, priceAll, pl, atAll, oAll, oTAll, oIndAll, tcAll, nSteps, datesStepRelAll, cbEffectiveDates, thetaF, Sigmaw, Sigmav, Sigmax, ax, Sigma0, Ef_)
        K.vStudentsT = true;

        K.nt = length(times);
        K.np = size(xpInit, 2);
        K.nc = nSteps;
        K.ns = K.np+K.nc;
        K.nu = sum(usedInstr);
        K.nx = K.np + K.nc + K.nu;
        K.nz = zeros(K.nt, 1);
        
        K.pl = pl;
        K.atAll = atAll;
        K.oAll = oAll;
        K.oTAll = oTAll;
        K.oIndAll = oIndAll;
        K.tcAll = tcAll;       
        K.datesStepRelAll = datesStepRelAll;
        
        K.x0 = zeros(K.nx,1);
        K.x = cell(K.nt,1);
        K.w = cell(K.nt,1);
        K.z = cell(K.nt,1);
        K.v = cell(K.nt,1);
        K.F = cell(K.nt,1);
        K.Iz = cell(K.nt,1);
        K.Cw = cell(K.nt,1);
        K.Cx = cell(K.nt,1);
        K.Cv = cell(K.nt,1);
        K.nuv = cell(K.nt,1);

        K.Ef = Ef_;
        K.EfInt = [zeros(1, size(K.Ef,2)) ; cumsum(K.Ef,1)]/365;

        % Parameters for Kalman filter

        K.thetaF = thetaF;
        K.Sigmaw = Sigmaw;
        K.Sigmav = Sigmav;
        K.Sigmax = Sigmax;
        K.ax = ax;
        K.Sigma0 = Sigma0;
        
        
%         K.thetaF = ones(K.nx,1);
%         K.Sigmaw = [ef*kx(1) ; ones(K.nc,1) ; ones(K.nu,1)*kx(2)];
%         K.Sigmav = ones(K.nu,1)*kx(4);
%         K.ax = zeros(K.nx,1);      % E[x0]
%         K.Sigma0 = eye(K.nx,1);   % Cov[x0] - denoted Sigmax by Carl & Axel - should be changed to Sigma0
%         K.Sigmax = [zeros(K.np,1) ; zeros(K.nc,1) ; ones(K.nu,1)*kx(3)];
        
        % Matrices for keeping track of parameters
        K.A = cell(K.nt,1); % ThetaF - Rows
        K.B = cell(K.nt,1); % ThetaF - Columns
        K.D = cell(K.nt,1); % Sigmaw
        K.G = cell(K.nt,1); % Sigmav
        
        
        te = find(cbEffectiveDates > floor(times(1)));
        te = te(1);
        
        for t=1:K.nt
          K.nz(t) = length(indInstrAll{t});
          K.x{t} = zeros(K.nx,1);
          K.w{t} = zeros(K.nx,1);
          K.z{t} = priceAll{t};
          K.v{t} = zeros(K.nz(t),1);
          if (floor(times(t)) >= cbEffectiveDates(te))
%             F = zeros(K.nx);
%             F(1:K.np, 1:K.np) = eye(K.np);
%             for i=1:K.nc-1
%               F(K.np+i, K.np+i+1) = 1;
%             end
%             F(K.ns+(1:K.nu), K.ns+(1:K.nu)) = eye(K.nu);
%             K.F{t} = sparse(F);
            te = te+1;
            K.A{t} = sparse(eye(K.nx));
            B = eye(K.nx);
            for i=(K.np+1):(K.np+K.nc)
              if (i < K.np+K.nc)
                B(i,i+1) = B(i,i);
              end
              B(i,i) = 0;
            end
            K.B{t} = sparse(B);
            
%             K.F{t} = sparse(eye(K.nx));
          else
%             K.F{t} = sparse(eye(K.nx));
            K.A{t} = sparse(eye(K.nx));
            K.B{t} = sparse(eye(K.nx));
          end
          K.Iz{t} = sparse(1:K.nz(t), indInstrAll{t}, ones(K.nz(t),1), K.nz(t), K.nu);

%           K.Cw{t} = [ef*kx(1) ; ones(K.nc,1) ; ones(K.nu,1)*kx(2)];
%           K.Cx{t} = [zeros(K.np,1) ; zeros(K.nc,1) ; ones(K.nu,1)*kx(3)];
          K.Cx{t} = K.Sigmax;
%           K.Cv{t} = ones(K.nz(t),1)*kx(4);
          if (K.vStudentsT)
            K.nuv{t} = ones(K.nz(t),1)*8;
          end
          K.x{t}(1:K.np) = xpInit(t,:)';
          
          F = K.A{t}*diag(K.thetaF)*K.B{t};
          K.F{t} = sparse(F);
%           if (norm(F-K.F{t}) > 1E-5)
%             error('Differ');
%           end
          
          K.G{t} = sparse(1:K.nz(t), indInstrAll{t}, ones(K.nz(t),1), K.nz(t), K.nu);
          Cv = K.G{t}*diag(K.Sigmav)*K.G{t}';
%           if (norm(diag(Cv)-K.Cv{t}) > 1E-5)
%             error('Differ');
%           end
          K.Cv{t} = diag(Cv);
          
          K.D{t} = sparse(eye(K.nx));
          Cw = K.D{t}*diag(K.Sigmaw)*K.D{t}';
%           if (norm(diag(Cw)-K.Cw{t}) > 1E-5)
%             error('Differ');
%           end
          K.Cw{t} = diag(Cw);
          
        end

      end % clsPriceStochastic constructor
      
      function [a, A, b, B, c, e, E, h, H, cxAll, cwAll, cvAll] = approximate(K, t, cxAll, cwAll, cvAll)
        % Second order Taylor approximation of normal pdf
        ind = (K.Cx{t}~=0);
      %   c = -sum(log(normpdf(K.w{t}, 0, sqrt(K.Cw{t}))))-sum(log(normpdf(K.x{t}(ind), 0, sqrt(K.Cx{t}(ind)))));
        c = -sum(-0.5*log(2*pi*K.Cw{t})-0.5*(K.w{t}./sqrt(K.Cw{t})).^2)-sum(-0.5*log(2*pi*K.Cx{t}(ind))-0.5*(K.x{t}(ind)./sqrt(K.Cx{t}(ind))).^2);
      cwAll = cwAll - sum(-0.5*log(2*pi*K.Cw{t})-0.5*(K.w{t}./sqrt(K.Cw{t})).^2);
      cxAll = cxAll - sum(-0.5*log(2*pi*K.Cx{t}(ind))-0.5*(K.x{t}(ind)./sqrt(K.Cx{t}(ind))).^2);

        if (K.vStudentsT) % Student's t-distribution
          vpdf = pdf('tLocationScale',K.v{t}, 0, sqrt(K.Cv{t}), K.nuv{t});
      %     vpdf2 = gamma((K.nuv{t}+1)/2)./(gamma(K.nuv{t}/2).*sqrt(K.nuv{t}*pi).*sqrt(K.Cv{t})).*(1 + K.v{t}.^2./(K.nuv{t}.*K.Cv{t})).^(-(K.nuv{t}+1)/2);
          c = c - sum(log(vpdf));
      cvAll = cvAll - sum(log(vpdf));
          k = gamma((K.nuv{t}+1)/2)./(gamma(K.nuv{t}/2).*sqrt(K.nuv{t}*pi).*sqrt(K.Cv{t}));
          dvpdf = -k.*(1 + K.v{t}.^2./(K.nuv{t}.*K.Cv{t})).^(-(K.nuv{t}+3)/2) .* (K.nuv{t}+1) .* K.v{t}./(K.nuv{t}.*K.Cv{t});
          d2vpdf = k.*(1 + K.v{t}.^2./(K.nuv{t}.*K.Cv{t})).^(-(K.nuv{t}+5)/2) .* (K.nuv{t}+3) .* (K.nuv{t}+1) .* K.v{t}.^2./(K.nuv{t}.^2.*K.Cv{t}.^2) ...
                    -k.*(1 + K.v{t}.^2./(K.nuv{t}.*K.Cv{t})).^(-(K.nuv{t}+3)/2) .* (K.nuv{t}+1)./(K.nuv{t}.*K.Cv{t});
          b = -dvpdf./vpdf;
          B = diag(dvpdf.^2./vpdf.^2 - d2vpdf./vpdf);

        else % Normal distribution
      %     c = c - sum(log(normpdf(K.v{t}, 0, sqrt(K.Cv{t}))));
          c = c - sum(-0.5*log(2*pi*K.Cv{t})-0.5*(K.v{t}./sqrt(K.Cv{t})).^2);
      cvAll = cvAll - sum(-0.5*log(2*pi*K.Cv{t})-0.5*(K.v{t}./sqrt(K.Cv{t})).^2);    
          b = K.v{t}./K.Cv{t};
          B = diag(1./K.Cv{t});
        end
        
        a = K.w{t}./K.Cw{t};
        A = diag(1./K.Cw{t});
        tmp = zeros(K.nx,1);
        tmp(ind) = K.x{t}(ind)./K.Cx{t}(ind);
        e = tmp;
        tmp(ind) = 1./K.Cx{t}(ind);
        E = diag(tmp);

        if (isinf(c))
          hej2 = 1
        end
        o = oCompute(K, t); % Instead of using K.oAll - allows updates of Ef
        [h, G] = taylorInstrPC(K.pl, K.atAll{t}, o, K.oIndAll{t}, K.tcAll{t}, K.x{t}(1:K.ns));
        H = [G K.Iz{t}];
      end

      function [f] = evalObj(K)
        % Evaluate objective function value
        f = 0; % Note that x0 term is missing
        K.EfInt = [zeros(1, size(K.Ef,2)) ; cumsum(K.Ef,1)]/365; % Ef can be changed             

        for t=1:K.nt         
          ind = (K.Cx{t}~=0);
%           f = f - sum(log(normpdf(K.w{t}, 0, sqrt(K.Cw{t}))))-sum(log(normpdf(K.x{t}(ind), 0, sqrt(K.Cx{t}(ind)))));
          f = f - sum(-0.5*log(2*pi*K.Cw{t})-0.5*(K.w{t}./sqrt(K.Cw{t})).^2)-sum(-0.5*log(2*pi*K.Cx{t}(ind))-0.5*(K.x{t}(ind)./sqrt(K.Cx{t}(ind))).^2);
          if (K.vStudentsT) % Student's t-distribution
            vpdf = pdf('tLocationScale',K.v{t}, 0, sqrt(K.Cv{t}), K.nuv{t});
            f = f - sum(log(vpdf));

          else % Normal distribution
%             f = f - sum(log(normpdf(K.v{t}, 0, sqrt(K.Cv{t}))));
            f = f - sum(-0.5*log(2*pi*K.Cv{t})-0.5*(K.v{t}./sqrt(K.Cv{t})).^2);
          end
        end
      end
      
      function [D] = dividends(obj, dm, dc)
        Nc = length(dm.cName);
        D = cell(Nc, 1);
      end

      function [iCurPrice] = priceCurrency(obj)
        iCurPrice = obj.iCurPrice;
      end
      
      function [f, KR] = logLikelihoodEps(K)
        % Evaluate log likelihood value
        f = 0; 
        KR.xPred = cell(K.nt,1);
        KR.PPred = cell(K.nt,1);
        KR.u = cell(K.nt,1);
        KR.g = cell(K.nt,1);
        KR.H = cell(K.nt,1);
        KR.eps = cell(K.nt,1);
        KR.S = cell(K.nt,1);
        KR.K = cell(K.nt,1);
        KR.xUpdate = cell(K.nt,1);
        KR.PUpdate = cell(K.nt,1);
        KR.fAll = zeros(K.nt,1);
        KR.fAll1 = zeros(K.nt,1);
        KR.fAll2 = zeros(K.nt,1);
        KR.fAll3 = zeros(K.nt,1);
        
        KR.tmp = cell(K.nt,1);
        K.EfInt = [zeros(1, size(K.Ef,2)) ; cumsum(K.Ef,1)]/365; % Ef can be changed             
        for t=1:K.nt
          F = sparse(K.A{t}*diag(K.thetaF)*K.B{t}); % Keep track of changes in thetaF
          Cv = K.G{t}*diag(K.Sigmav)*K.G{t}'; % Keep track of changes in Sigmav
          Cw = K.D{t}*diag(K.Sigmaw)*K.D{t}'; % Keep track of changes in Sigmaw

          if (t==1)
            KR.xPred{t} = F * K.ax;
            KR.PPred{t} = F * diag(K.Sigma0) * F' + Cw;
          else
            KR.xPred{t} = F * KR.xUpdate{t-1};
            KR.PPred{t} = F * KR.PUpdate{t-1} * F' + Cw;
          end
          o = oCompute(K, t); % Instead of using K.oAll - allows updates of Ef
          [g, G] = taylorInstrPC(K.pl, K.atAll{t}, o, K.oIndAll{t}, K.tcAll{t}, KR.xPred{t}(1:K.ns));
          KR.u{t} = g - G*KR.xPred{t}(1:K.ns);
          KR.g{t} = g;
          KR.H{t} = [G K.Iz{t}];
          KR.eps{t} = K.z{t} - KR.H{t}*KR.xPred{t} - KR.u{t};
          PH = KR.PPred{t} * KR.H{t}';
          KR.S{t} = KR.H{t} * PH + Cv;
%           KR.K{t} = PH/KR.S{t};
          KR.K{t} = PH*inv(KR.S{t});
          KR.xUpdate{t} = KR.xPred{t} + KR.K{t}*KR.eps{t};
%           KR.PUpdate{t} = KR.PPred{t} - KR.K{t}*PH';
          KR.PUpdate{t} = KR.PPred{t} - KR.K{t}*KR.H{t}*KR.PPred{t};
          f = f - 0.5*(K.nz(t)*log(2*pi) + log(det(KR.S{t})) + KR.eps{t}' * (KR.S{t} \ KR.eps{t}));
          KR.fAll(t) = - 0.5*(K.nz(t)*log(2*pi) + log(det(KR.S{t})) + KR.eps{t}' * (KR.S{t} \ KR.eps{t}));
          KR.fAll1(t) = - 0.5*K.nz(t)*log(2*pi);
          KR.fAll2(t) = - 0.5*log(det(KR.S{t}));
          KR.fAll3(t) = - 0.5* KR.eps{t}' * (KR.S{t} \ KR.eps{t});
KR.tmp{t} = g;
% KR.tmp{t} = KR.S{t};
          if (t>1)
% KR.tmp{t} = K.B{t} * KR.PUpdate{t-1} * K.B{t}';
          end
          if (false) % Compute implied v from xUpdate
            [g, G] = taylorInstrPC(K.pl, K.atAll{t}, o, K.oIndAll{t}, K.tcAll{t}, KR.xUpdate{t}(1:K.ns));
            KR.tmp{t} = K.z{t} - g - K.Iz{t}*KR.xUpdate{t}(K.ns+1:end);
          end
        end
      end

      function [] = storeParamVec(K, paramVec)
        thetaFInd = (1:length(K.thetaF));
        SigmavInd = (1:length(K.Sigmav))+thetaFInd(end);
        SigmawInd = (1:length(K.Sigmaw))+SigmavInd(end);
        axInd = (1:length(K.ax))+SigmawInd(end);
        Sigma0Ind = (1:length(K.Sigma0))+axInd(end);
        EfVecInd = (1:length(K.Ef(:)))+Sigma0Ind(end);

        K.thetaF = paramVec(thetaFInd);
        K.Sigmaw = paramVec(SigmawInd);
        K.Sigmav = paramVec(SigmavInd);
        K.ax = paramVec(axInd);
        K.Sigma0 = paramVec(Sigma0Ind);
        K.Ef(:) = paramVec(EfVecInd);
      end      
      
      function [xPred, xUpdate, eps, z, g] = analyzeDetermineStates(K, KI, paramVec)
        K.storeParamVec(paramVec);
        [f, KR] = K.logLikelihoodEps();
        
        xPred = zeros(K.nt, K.nx);
        xUpdate = zeros(K.nt, K.nx);
        eps = zeros(K.nt, K.nu);
        z = nan(K.nt, K.nu);
        g = nan(K.nt, K.nu);
        for t=1:K.nt
          xPred(t,:) = KR.xPred{t};
          xUpdate(t,:) = KR.xUpdate{t};
          eps(t,KI.indInstrAll{t}) = KR.eps{t};
          z(t,KI.indInstrAll{t}) = K.z{t};
          g(t,KI.indInstrAll{t}) = KR.g{t};
        end
      end

      function [xp, xc, xu] = analyzeSplitState(K, x)
        xp = x(:,1:K.np);
        xc = x(:,K.np + (1:K.nc));
        xu = x(:,K.ns+1:end);
      end
      
      function [eps2] = analyze(K, KI, paramVec)
        K.storeParamVec(paramVec);
        [f, KR] = K.logLikelihoodEps();
        fprintf('Log-likelihood = %d\n', f);

        [xPred, xUpdate, eps, z, g] = K.analyzeDetermineStates(KI, paramVec);

        [xp, xc, xu] = K.analyzeSplitState(xPred);

        figure(11);
        plot(xp(:,1:3));
        hold on;
        plot(xp(:,4:6), ':');        
        hold off;
        title('x PC');

        figure(12);
        plot(xc(:,1:3));
        title('x step');

        figure(13);
        plot(xu);
        title('x unsystematic');

        [xUp, xUc, xUu] = K.analyzeSplitState(xUpdate);

        figure(14);
        plot(xUp(:,1:3));
        hold on;
        plot(xUp(:,4:6), ':');        
        hold off;
        title('xU PC');

        figure(15);
        plot(xUc(:,1:3));
        title('xU step');

        figure(16);
        plot(xUu);
        title('xU unsystematic');


        figure(17);
        plot(eps);
        title('epsilon');
        
        eps2 = sum(eps.^2,1);
      
        figure(18);
        plot(KR.fAll2);
        title('fAll2');
        
        figure(19);
        plot(KR.fAll3);
        title('fAll3');

        fprintf('%15s %13s %13s\n', 'Risk factor', 'Sigmaw', 'var(x)');
        for i=1:K.nx
          if (i <= K.np)
            xStr = strcat('PC', num2str(i));
            varx = mean(xp(:,i).^2);
          elseif (i <= K.np+K.nc)
            ii = i - K.np;
            varx = mean(xc(:,ii).^2);
            xStr = strcat('Cut/Hike', num2str(ii));
          else
            ii = i - K.ns;
            varx = mean(xu(:,ii).^2);
            xStr = KI.ricAll{ii};
          end
          fprintf('%15s %13d %13d\n', xStr, K.Sigmaw(i), varx);
        end
        
        fprintf('%15s %13s %13s\n', 'RIC', 'Sigmav', 'var(xu)');
        for i=1:K.nu
          fprintf('%15s %13d %13d\n', KI.ricAll{i}, K.Sigmav(i), mean(xu(:,i).^2));
        end

        i = 3;
        figure(20);
        plot([z(:,i) g(:,i)])
        title(KI.ricAll{i});
      end
      
      function [] = plotSolution(K, KI, paramVec)
        K.storeParamVec(paramVec);

        figure(1);
        X = 1:length(K.thetaF);
        plot(X, K.thetaF, '-', X, K.thetaF, ':');
        title('thetaF');

        figure(2);
        X = 1:length(K.Sigmaw);
        semilogy(X, K.Sigmaw, '-', X, K.Sigmaw, ':');
        title('Sigmaw');

        figure(3);
        X = 1:length(K.Sigmav);
        semilogy(X, K.Sigmav, '-', X, K.Sigmav, ':');
        title('Sigmav');

        figure(4);
        X = 1:length(K.ax);
        plot(X, K.ax, '-', X, K.ax, ':');
        title('ax');

        figure(5);
        X = 1:length(K.Sigma0);
        semilogy(X, K.Sigma0, '-', X, K.Sigma0, ':');
        title('Sigma0');

        figure(6);
        T = (0:size(K.Ef,1)-1)/365;
        subplot(2,1,1);
        plot(T, K.Ef(:,1:3));
        title('Ef(1-3)');
        subplot(2,1,2);
        plot(T, K.Ef(:,4:6));
        title('Ef(4-6)');
      end
      
      function [] = plotCurves(K, KI, paramVec)
        K.storeParamVec(paramVec);
        [f, KR] = K.logLikelihoodEps();

        [xPred, xUpdate, eps, z, g] = K.analyzeDetermineStates(KI, paramVec);

        [xp, xc, xu] = K.analyzeSplitState(xPred);

        for k=4764:K.nt
          datesStepRel = KI.datesStepRelAll{k};
          Es = zeros(size(K.Ef,1), K.nc);
          for i=1:K.nc
            Es(datesStepRel(i)+1:end, i) = 1;
          end
          E = [K.Ef Es];
          
          plot((0:size(E,1)-1)/365, E*[xp(k,:) xc(k,:)]');
          title(datestr(floor(KI.times(k))));
          pause(0.01);

        end
      end
      
   end
   
   
   methods (Access = 'private') % Access by class members only
      function [o] = oCompute(K, t)
        datesStepRel = K.datesStepRelAll{t};
        o = [-K.EfInt(K.oTAll{t},:) zeros(size(K.oTAll{t},1), K.nc)];
        for i=1:K.nc
          o(:,K.np+i) = -max(K.oTAll{t}-(datesStepRel(i)+1), 0)/365;
        end
      end
           
   end
end % classdef





