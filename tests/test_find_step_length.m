
% test_full_line_search_pipeline.m
% Testar hela pipelinen:
% poly_fit_ml -> mult_q_poly_fits -> select_poly_degree -> find_step_length
%
% Testar funktioner som:
% 1. Ökar initialt och sedan når ett maximum.
% 2. Är avtagande redan från start.
%
% För varje testfall visas:
% - Mockdata
% - Anpassad ML-kurva
% - Vald steglängd och motsvarande funktionsvärde

clear;
clc;
close all;

%% Lägg till sökvägar till funktionerna

projectRoot = fileparts(fileparts(mfilename('fullpath')));

addpath(fullfile(projectRoot, 'scripts', 'line_search'));

%% Inställningar

rng(42);

m = 8;              % Antal datapunkter
qMax = 5;             % Högsta polynomgrad som testas

noiseF = 0.1;        % Brus i funktionsvärden
noiseS = 0.1;        % Brus i derivatan

alpha = linspace(0, 1.5, m)';

% Parametrar till find_step_length
alphaMax = 1.5;

% Vi väljer ett initialt steg som vi vet ger förbättring.
% På så sätt testar vi främst hur find_step_length hittar maximum.
alphaPrev = alpha(2);
isBetter = true;

%% Definiera testfunktioner

% Varje funktion har en tillhörande derivata.
% Derivatan motsvarar S i poly_fit_ml.

cases(1).name = 'Andragradspolynom';
cases(1).F = @(a) 1 + 4*a - 2*a.^2;
cases(1).S = @(a) 4 - 4*a;

cases(2).name = 'Tredjegradspolynom';
cases(2).F = @(a) 1 + 3*a + 0.5*a.^2 - 1.5*a.^3;
cases(2).S = @(a) 3 + a - 4.5*a.^2;

cases(3).name = 'Fjärdegradspolynom';
cases(3).F = @(a) 2 + 5*a - a.^2 - 0.5*a.^3 - 0.2*a.^4;
cases(3).S = @(a) 5 - 2*a - 1.5*a.^2 - 0.8*a.^3;

% Dessa funktioner är avtagande redan vid alpha = 0.
cases(4).name = 'Avtagande linjär funktion';
cases(4).F = @(a) 5 - 2*a;
cases(4).S = @(a) -2 + 0*a;

cases(5).name = 'Avtagande andragradspolynom';
cases(5).F = @(a) 4 - 2*a - 0.5*a.^2;
cases(5).S = @(a) -2 - a;

%% Kör testerna

for test = 1:length(cases)

    fprintf('\n');
    fprintf('============================================\n');
    fprintf('TEST %d: %s\n', test, cases(test).name);
    fprintf('============================================\n');

    try

        %% 1. Generera mockdata

        F_true = cases(test).F(alpha);
        S_true = cases(test).S(alpha);

        F = F_true + noiseF*randn(m,1);
        S = S_true + noiseS*randn(m,1);

        fprintf('Data genererad.\n');

        %% 2. Kör mult_q_poly_fits
        %
        % Denna funktion kör poly_fit_ml för alla grader
        % från 1 till qMax och beräknar residualdiagnostik.

        fprintf('\nKör mult_q_poly_fits...\n');

        fits = mult_q_poly_fits(alpha, F, S, qMax);

        %% 3. Välj polynomgrad

        fprintf('\nKör select_poly_degree...\n');

        [q, fit, info] = select_poly_degree(fits, 'threshold');

        % Om tröskelmetoden inte hittar någon godkänd grad
        % testar vi den poängbaserade metoden.

        if isempty(q)
            fprintf(['Tröskelmetoden hittade ingen godkänd grad. ' ...
                     'Testar score-metoden...\n']);

            [q, fit, info] = select_poly_degree(fits, 'score');
        end

        if isempty(q)
            error('Ingen användbar polynomgrad kunde väljas.');
        end

        fprintf('Vald polynomgrad: %d\n', q);

        %% 4. Kör find_step_length

        fprintf('\nKör find_step_length...\n');

        [alphaNew, alphaMaxNew] = find_step_length( ...
            fit, alphaPrev, alphaMax, isBetter);

        % Beräkna funktionsvärdet från den anpassade ML-kurvan.
        FNew = polyval(flipud(fit.y(:)), alphaNew);

        fprintf('Vald steglängd: %.4f\n', alphaNew);
        fprintf('Funktionsvärde vid vald steglängd: %.4f\n', FNew);
        fprintf('Uppdaterad alphaMax: %.4f\n', alphaMaxNew);

        %% 5. Kontrollera om funktionen ökar initialt

        initialDerivative = polyval( ...
            flipud((1:q)' .* fit.y(2:end)), 0);

        if initialDerivative <= 0
            fprintf('\nVARNING: Den anpassade kurvan är inte ');
            fprintf('ökande vid alpha = 0.\n');
            fprintf(['find_step_length har ändå returnerat en ' ...
                     'steglängd enligt sin nuvarande logik.\n']);
        else
            fprintf('\nDen anpassade kurvan ökar vid alpha = 0.\n');
        end

        %% 6. Visualisera resultatet

        alphaGrid = linspace(0, alphaMax, 300)';

        % Anpassad ML-kurva
        F_fit = polyval(flipud(fit.y(:)), alphaGrid);

        % ML-kurvans värde vid vald steglängd
        F_marker = polyval(flipud(fit.y(:)), alphaNew);

        figure('Name', cases(test).name);

        hold on;
        grid on;

        % Mockdata
        scatter(alpha, F, 35, 'b', 'filled', ...
            'DisplayName', 'Mockdata');

        % Anpassad ML-kurva
        plot(alphaGrid, F_fit, 'r-', 'LineWidth', 2, ...
            'DisplayName', sprintf('ML-anpassning (grad %d)', q));

        % Vald steglängd
        plot(alphaNew, F_marker, 'ko', ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', 'g', ...
            'DisplayName', 'Vald steglängd');

        % Vertikal linje vid vald steglängd
        xline(alphaNew, 'k--', ...
            'HandleVisibility', 'off');

        % Markera funktionsvärdet
        yline(F_marker, 'k:', ...
            'HandleVisibility', 'off');

        % Text vid markeringen
        text(alphaNew, F_marker, ...
            sprintf('  alpha = %.3f\\n  F = %.3f', ...
            alphaNew, F_marker), ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9);

        xlabel('\alpha');
        ylabel('F(\alpha)');

        title(sprintf('%s - vald grad %d', ...
            cases(test).name, q));

        legend('Location', 'best');

        hold off;

        fprintf('\nTest %d slutfört.\n', test);

    catch ME

        %% Felhantering

        fprintf('\nTEST %d MISSLYCKADES!\n', test);
        fprintf('Felmeddelande: %s\n', ME.message);

        % Fortsätt med nästa testfall.
        continue;

    end

end

fprintf('\n============================================\n');
fprintf('ALLA TESTFALL HAR BEHANDLATS.\n');
fprintf('============================================\n');
