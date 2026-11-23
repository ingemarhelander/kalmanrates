% replayDegreeSelection.m
% Spelar upp gradvalet från en riktig körning: kör select_poly_degree med
% och utan hoppet q -> q+2 (opts.lookAhead) på exakt samma punkter som
% linjesökningen använde, och plottar skillnaden. Inga Kalman-körningar.
%
% Kör först kalmanOptimizeMain med ls_saveFits 1 i configfilen. Den sparar
% data/degreeLog_<tid>.mat. Utan logFile nedan används den senaste loggen.
%
% Figurer:
%   1. Vald grad för varje gradval (gammalt mot nytt) och hur långt
%      polynomets maximum flyttas.
%   2. LR för varje test mot gränserna, mot antal punkter m.
%   3. Ett enskilt gradval (showFit, standard: det där maximum flyttas
%      mest): punkterna med polynomen för gammalt och nytt gradval.

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'line_search'));

logFile = '';   % tom: senaste data/degreeLog_*.mat
showFit = [];   % löpnummer för figur 3, tom: där maximum flyttas mest

%% Ladda loggen
if isempty(logFile)
    files = dir(dataPath('degreeLog_*.mat'));
    if isempty(files)
        error('kalmanrates:NoDegreeLog', ...
              'Ingen data/degreeLog_*.mat. Sätt ls_saveFits 1 och kör kalmanOptimizeMain först.');
    end
    [~, newest] = max([files.datenum]);
    logFile = fullfile(files(newest).folder, files(newest).name);
end
L = load(logFile);
config = readConfig(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));
h = 1e-3;
if isfield(config, 'ls_h')
    h = config.ls_h;
end
fprintf('Logg: %s (modell %s, %d linjesökningar)\n', logFile, L.modelFile, numel(L.searches));

%% Kör gammalt och nytt gradval på samma punkter
rows = struct('search', {}, 'iter', {}, 'denoise', {}, 'k', {}, 'm', {}, ...
              'alphaMax', {}, 'qOld', {}, 'qNew', {}, 'tHatOld', {}, 'tHatNew', {}, ...
              'fitOld', {}, 'fitNew', {}, 'infoNew', {}, 'data', {});
for s = 1:numel(L.searches)
    S = L.searches(s);
    for j = 1:numel(S.fits)
        d = S.fits{j};
        % evalc: poly_fit_ml skriver ut när den inte konvergerar
        [~, qOld, fitOld] = evalc('select_poly_degree(d.t, d.F, d.St, struct(''lookAhead'', false))');
        [~, qNew, fitNew, infoNew] = evalc('select_poly_degree(d.t, d.F, d.St, struct(''lookAhead'', true))');
        rows(end+1) = struct('search', s, 'iter', S.iter, 'denoise', S.denoise, ...
                             'k', d.k, 'm', numel(d.t), 'alphaMax', d.alphaMax, ...
                             'qOld', qOld, 'qNew', qNew, ...
                             'tHatOld', find_step_length(fitOld, h, 1), ...
                             'tHatNew', find_step_length(fitNew, h, 1), ...
                             'fitOld', fitOld, 'fitNew', fitNew, 'infoNew', infoNew, ...
                             'data', d); %#ok<SAGROW>
    end
end
if isempty(rows)
    error('kalmanrates:EmptyDegreeLog', '%s innehåller inga gradval.', logFile);
end

qOld = [rows.qOld];
qNew = [rows.qNew];
shift = abs([rows.tHatNew] - [rows.tHatOld]);   % flytt av maximum, relativt alphaMax
changed = qOld ~= qNew;
fprintf('%d gradval, grad ändrad i %d (%.0f %%), maximum flyttat > 1 %% av alphaMax i %d\n', ...
        numel(rows), sum(changed), 100*mean(changed), sum(shift > 0.01));

%% Figur 1: vald grad och flytt av maximum
n = 1:numel(rows);
newSearch = find(diff([rows.search]) ~= 0) + 0.5;   % gränser mellan linjesökningar
figure('Name', 'Gradval per anpassning');
subplot(2, 1, 1);
hold on;
stairs([n, n(end) + 1] - 0.5, [qOld, qOld(end)], 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);
plot(n, qNew, 'o', 'MarkerSize', 4);
plot(n(changed), qNew(changed), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
xline(newSearch, ':', 'Color', [0.7 0.7 0.7]);
ylabel('grad q');
ylim([1.5 6.5]);
legend({'utan q \rightarrow q+2', 'med q \rightarrow q+2', 'ändrad'}, 'Location', 'best');
title('Vald grad (prickade linjer: ny linjesökning)');
subplot(2, 1, 2);
semilogy(n, max(shift, 1e-6), '.-');
hold on;
xline(newSearch, ':', 'Color', [0.7 0.7 0.7]);
yline(0.01, '--');
xlabel('gradval nr (alla linjesökningar i följd)');
ylabel('|\Deltat_{max}| (relativt alphaMax)');
title('Hur långt polynomets maximum flyttas (golv 1e-6)');

%% Figur 2: LR mot gränserna
LR = [];  m = [];  jump = [];  raised = [];
for r = 1:numel(rows)
    I = rows(r).infoNew;
    LR = [LR, I.LR]; %#ok<AGROW>
    m = [m, repmat(I.m, 1, numel(I.LR))]; %#ok<AGROW>
    jump = [jump, I.qTo - I.qFrom]; %#ok<AGROW>
    raised = [raised, I.LR > I.critTest]; %#ok<AGROW>
end
figure('Name', 'LR-test');
hold on;
one = jump == 1;
scatter(m(one) + 0.1*randn(1, sum(one)), max(LR(one), 1e-2), 15, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(m(~one) + 0.1*randn(1, sum(~one)), max(LR(~one), 1e-2), 30, 'd', 'filled');
yline(rows(1).infoNew.crit, '--', 'crit (q\rightarrowq+1)');
yline(rows(1).infoNew.crit2, '-.', 'crit2 (q\rightarrowq+2)');
set(gca, 'YScale', 'log');
xlabel('antal punkter m');
ylabel('LR (golv 1e-2)');
legend({'q \rightarrow q+1', 'q \rightarrow q+2'}, 'Location', 'best');
title(sprintf('LR per test: %d av %d höjer graden', sum(raised), numel(LR)));

%% Figur 3: ett enskilt gradval
if isempty(showFit)
    [~, showFit] = max(shift);
end
R = rows(showFit);
d = R.data;
alphaGrid = linspace(0, max(d.t), 400)*d.alphaMax;
tGrid = alphaGrid/d.alphaMax;
F0 = d.F(d.t == 0);
if isempty(F0)
    F0 = 0;
end
F0 = F0(1);
figure('Name', sprintf('Gradval nr %d', showFit));
subplot(2, 1, 1);
plot(d.t*d.alphaMax, d.F - F0, 'ko', 'MarkerFaceColor', 'k');
hold on;
plot(alphaGrid, poly_value(R.fitOld, tGrid) - F0, '--', 'LineWidth', 1.5);
plot(alphaGrid, poly_value(R.fitNew, tGrid) - F0, '-', 'LineWidth', 1.5);
xline(R.tHatOld*d.alphaMax, '--');
xline(R.tHatNew*d.alphaMax, '-');
ylabel('F - F_0');
legend({'mätt', sprintf('grad %d (utan hopp)', R.qOld), sprintf('grad %d (med hopp)', R.qNew)}, ...
       'Location', 'best');
title(sprintf('Iteration %d%s, punkt %d, m = %d, alphaMax = %.3g', R.iter, ...
              repmat(' (denoise)', 1, R.denoise), R.k, R.m, d.alphaMax));
subplot(2, 1, 2);
% Lutningen i alpha: dF/dalpha = (dF/dt)/alphaMax
plot(d.t*d.alphaMax, d.St/d.alphaMax, 'ko', 'MarkerFaceColor', 'k');
hold on;
plot(alphaGrid, poly_slope(R.fitOld, tGrid)/d.alphaMax, '--', 'LineWidth', 1.5);
plot(alphaGrid, poly_slope(R.fitNew, tGrid)/d.alphaMax, '-', 'LineWidth', 1.5);
yline(0, ':');
xlabel('alpha');
ylabel('dF/d\alpha');
