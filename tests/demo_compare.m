% demo_compare - Jämför gradval med lag-1 och LR-test för få punkter.
%
% För varje funktion och antal punkter m slumpas 50 dataset, och båda
% metoderna väljer grad. Tabellen visar andelen dataset där rätt grad
% valdes, eller medelvärdet av vald grad för funktioner som inte är polynom.
%
% Kör från projektroten med:  demo_compare

clear;
clc;

addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'scripts', 'line_search'));

sF = 0.003;   % brus på F
sS = 0.015;   % brus på S
ms = [4 5 6 8 10];
nTrials = 50;

% Namn, sann F, sann S, sann grad (NaN = inget polynom)
funs = {
    'Kvadratisk',    @(a) 1 + 3*a - 5*a.^2,    @(a) 3 - 10*a,             2
    'Kubisk',        @(a) 1 + 2*a - 6*a.^3,    @(a) 2 - 18*a.^2,          3
    'Fjärdegrads',   @(a) 1 - 4*a.^3 + 5*a.^4, @(a) -12*a.^2 + 20*a.^3,   4
    '1 - exp(-4a)',  @(a) 1 - exp(-4*a),       @(a) 4*exp(-4*a),          NaN
};
selectors = {'lag1', 'lr'};

fprintf('%-14s %4s   %-10s %-10s\n', 'Funktion', 'm', 'lag-1', 'LR');
for i = 1:size(funs, 1)
    for m = ms
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
        if isnan(funs{i, 4})
            fprintf('%-14s %4d   q = %-6.1f q = %-6.1f\n', funs{i, 1}, m, mean(qq(:, 1)), mean(qq(:, 2)));
        else
            fprintf('%-14s %4d   %5.0f %%    %5.0f %%\n', funs{i, 1}, m, 100*mean(qq == funs{i, 4}));
        end
    end
    fprintf('\n');
end
