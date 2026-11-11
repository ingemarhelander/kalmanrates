function fit = poly_fit_ml(alpha, F, S, q, doPlot)
% poly_fit_ml - Total ML-skattning av ett polynom med både
% funktionsvärden och riktningsderivator.
%
% Skapar polynom givet gradtal q, steglängd alpha, funktionsvärden F och riktningsderivator S.
%
% Implementerar den totala ML-skattningen enligt ekvation
% (34), (38) och (39) i ML-dokumentet.
%
% In:
%   alpha  Steglängder för de olika observationerna (m x 1)
%   F      Funktionsvärden (m x 1)
%   S      Lutningar längs linjen, g'*d (m x 1)
%   q      Polynomgrad
%  doPlot  (valfritt) om true, rita residualer och konvergens (TA BORT EFTER TESTNING)
%
% Ut:
%   fit.y       Polynomkoefficienter (q+1 x 1)
%   fit.rF      Residualer för funktionsvärden (m x 1)
%   fit.rG      Residualer för lutningar (m x 1)
%   fit.wF      Vikt för F = 1/sigmaF^2
%   fit.wG      Vikt för S = 1/sigmaG^2
%   fit.sigmaF  Skattad standardavvikelse för F
%   fit.sigmaG  Skattad standardavvikelse för S
%
% ML-metod:
%   1. Initiera vikterna wF = wG = 1.
%   2. Skatta y genom gemensam viktad minsta kvadrat.
%   3. Beräkna residualerna.
%   4. Skatta brusvarianserna sigmaF^2 och sigmaG^2.
%   5. Uppdatera vikterna.
%   6. Upprepa från steg 2 tills konvergens.

    % Ta bort när poly_fit_ml är testad och fungerar. Detta gör att plottar genereras
    if nargin < 5
        doPlot = true;
    end

    % Säkerställ kolumnvektorer
    alpha = alpha(:);
    F = F(:);
    S = S(:);

    m = length(alpha);

    % Check input dimensions
    if length(F) ~= m || length(S) ~= m
        error('All input vectors must have the same number of observations.');
    end

    if q < 1 || q ~= floor(q)
        error('q must be a positive integer.');
    end

    % Polynommatris och dess derivata
    a = bsxfun(@power, alpha, 0:q);

    b = [zeros(m,1), ...
         bsxfun(@times, bsxfun(@power, alpha, 0:q-1), 1:q)];


    % Check that the polynomial coefficients are identifiable (KAN TAS BORT EFTER TESTNING)
    if rank([a; b]) < q+1
        error(['The design matrix is rank deficient. ', ...
               'The selected step lengths cannot identify all coefficients.']);
    end
    
    % Initialisering enligt dokumentet
    wF = 1;
    wG = 1;

    % Iterationsparametrar
    maxIter = 1000;
    tol = 1e-8;
    converged = false;

    floorF = eps(max(sqrt(mean(F.^2)), 1))^2;
    floorG = eps(max(sqrt(mean(S.^2)), 1))^2;

    for k = 1:maxIter

        % Weighted ML estimate of polynomial coefficients
        A = [sqrt(wF)*a;
             sqrt(wG)*b];

        rhs = [sqrt(wF)*F;
               sqrt(wG)*S];

        y = A \ rhs;

        % Residualer för respektive mätserie
        rF = F - a*y;
        rG = S - b*y;

        % ML-skattning av brusvarianserna.
        % Motsvarar ekvation (38) och (39).
        sigmaF2 = max(mean(rF.^2), floorF);
        sigmaG2 = max(mean(rG.^2), floorG);

        % Uppdatera vikterna
        wF_new = 1/sigmaF2;
        wG_new = 1/sigmaG2;

        % Kontrollera konvergens
        change = max(abs([wF_new/wF - 1, ...
                          wG_new/wG - 1]));

        wF = wF_new;
        wG = wG_new;

        % Kontrollera konvergens
        if change < tol
            converged = true;
            break
        end
    end

    if converged
        disp('ML iteration converged')
    else
        disp('Maximum number of iterations reached without convergence')
    end

    % Beräkna om uppskattningarna med de slutliga vikterna
    A = [sqrt(wF)*a;
         sqrt(wG)*b];

    rhs = [sqrt(wF)*F;
           sqrt(wG)*S];

    y = A \ rhs;

    % Slutliga residualer och varianser
    rF = F - a*y;
    rG = S - b*y;

    sigmaF2 = max(mean(rF.^2), floorF);
    sigmaG2 = max(mean(rG.^2), floorG);

    % Resultat
    fit.y = y;

    fit.rF = rF;
    fit.rG = rG;

    fit.sigmaF = sqrt(sigmaF2);
    fit.sigmaG = sqrt(sigmaG2);

    fit.wF = 1/sigmaF2;
    fit.wG = 1/sigmaG2;

    % Information om iterationen
    fit.iterations = k;
    fit.converged = (k < maxIter);

    % Visualize measurements and fitted ML polynomial
    if doPlot

        alphaPlot = linspace(min(alpha), max(alpha), 300)';

        % Evaluate fitted polynomial
        Fplot = polyval(flipud(y), alphaPlot);

        % Evaluate derivative of fitted polynomial
        dy = (1:q)' .* y(2:end);
        Splot = polyval(flipud(dy), alphaPlot);

        figure('Name', 'ML polynomial fit');

        % Function values
        subplot(2,1,1);
        plot(alpha, F, 'ko', 'MarkerFaceColor', 'b');
        hold on;
        plot(alphaPlot, Fplot, 'r-', 'LineWidth', 1.5);
        xlabel('Step length \alpha');
        ylabel('Function value F');
        title(['ML polynomial fit, q = ', num2str(q)]);
        legend('Measurements', 'ML polynomial', 'Location', 'best');
        grid on;
        hold off;

        % Derivative values
        subplot(2,1,2);
        plot(alpha, S, 'ko', 'MarkerFaceColor', 'b');
        hold on;
        plot(alphaPlot, Splot, 'r-', 'LineWidth', 1.5);
        xlabel('Step length \alpha');
        ylabel('Derivative S');
        title('Derivative fit');
        legend('Measurements', 'ML derivative', 'Location', 'best');
        grid on;
        hold off;

    end

end