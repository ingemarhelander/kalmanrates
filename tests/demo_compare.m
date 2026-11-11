% demo_compare - Jämför gradval med lag-1 och LR-test för få punkter.
%
% För varje funktion och antal punkter m = 2..10 slumpas 50 dataset, och
% båda metoderna väljer grad. En figur per funktion visar medelvärdet av
% vald grad mot m, tillsammans med den sanna graden. En tabell i
% kommandofönstret visar andelen dataset där rätt grad valdes.
%
% Kör från projektroten med:  demo_compare

clear;
clc;
close all;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

sF = 0.003;   % brus på F
sS = 0.015;   % brus på S
ms = 2:10;
nTrials = 50;

% Namn, sann F, sann S, sann grad
funs = {
    'Kvadratisk',    @(a) 1 + 3*a - 5*a.^2,    @(a) 3 - 10*a,             2
    'Kubisk',        @(a) 1 + 2*a - 6*a.^3,    @(a) 2 - 18*a.^2,          3
    'Fjärdegrads',   @(a) 1 - 4*a.^3 + 5*a.^4, @(a) -12*a.^2 + 20*a.^3,   4
};
selectors = {'lag1', 'lr'};

% Färger som syns på både ljus och mörk bakgrund
orange = [1.0 0.6 0.0];
blue   = [0.3 0.6 1.0];
gray   = [0.6 0.6 0.6];

fprintf('Andel dataset där rätt grad valdes\n');
fprintf('%-14s %4s   %-8s %-8s\n', 'Funktion', 'm', 'lag-1', 'LR');
for i = 1:size(funs, 1)
    meanQ = zeros(numel(ms), 2);
    for j = 1:numel(ms)
        m = ms(j);
        alpha = linspace(0, 1, m)';
        qq = zeros(nTrials, 2);
        for s = 1:nTrials
            rng(s);
            F = funs{i, 2}(alpha) + sF*randn(m, 1);
            S = funs{i, 3}(alpha) + sS*randn(m, 1);
            for k = 1:2
                qk = [];
                evalc('qk = ml_poly_main(alpha, F, S, [], selectors{k});');   % utan utskrifter
                qq(s, k) = qk;
            end
        end
        meanQ(j, :) = mean(qq, 1);
        fprintf('%-14s %4d   %5.0f %%  %5.0f %%\n', funs{i, 1}, m, 100*mean(qq == funs{i, 4}));
    end
    fprintf('\n');

    % Vald grad mot antal punkter
    figure('Name', funs{i, 1});
    hold on;
    plot(ms, meanQ(:, 1), 'o-', 'Color', orange, 'MarkerFaceColor', orange, 'LineWidth', 2);
    plot(ms, meanQ(:, 2), 's-', 'Color', blue, 'MarkerFaceColor', blue, 'LineWidth', 2);
    yline(funs{i, 4}, '--', 'Color', gray, 'LineWidth', 1.5);
    hold off;
    grid on;
    xlim([ms(1) - 0.5, ms(end) + 0.5]);
    ylim([1.5 5.5]);
    xticks(ms);
    yticks(2:5);
    xlabel('Antal punkter m');
    ylabel(sprintf('Vald grad (medel av %d dataset)', nTrials));
    legend({'lag-1', 'LR', 'Sann grad'}, 'Location', 'northwest');
    title(sprintf('%s: sann grad %d', funs{i, 1}, funs{i, 4}));
end
