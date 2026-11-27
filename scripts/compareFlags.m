% compareFlags.m
% Are the quotes flagged with the START parameters (what prefilter in
% kalmanOptimizeMain does) the same as those flagged with OPTIMIZED
% parameters? If they mostly agree, the short optimization before the
% cleaning is not needed.
%
% Uses the model of kalmanOptimizeMain.config and flags, from the first
% optimization day to the last day of the model:
%   A  start parameters of the model (buildKalmanModel), robust scale
%   B  optimized parameters (data/kalmanOptimizeResult.mat), robust scale
%   C  optimized parameters, absolute scale (as cleanMarketData)
% B and C need a kalmanOptimizeResult.mat made on the same model; without
% it only A is shown.
%
% Settings below can be set before the script is run, e.g.
%   threshold = 5; compareFlags

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));

%% Settings
if ~exist('threshold', 'var'), threshold = 4; end
if ~exist('burnIn', 'var'), burnIn = 20; end
if ~exist('maxFlagFraction', 'var'), maxFlagFraction = 0.2; end
if ~exist('nExamplePlots', 'var'), nExamplePlots = 4; end

%% Model and parameters
cmpConfig = readConfig(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));
cmpPeriods = experimentPeriods(cmpConfig);
[K, times, ricAll, KI] = loadKalmanModel(cmpPeriods.modelFile);
optimizeIdx = selectDays(times, cmpPeriods.optimizeStart, cmpPeriods.optimizeEnd);
days = optimizeIdx(1):K.nt;
xStart = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];

opts = struct('threshold', threshold, 'burnIn', burnIn, 'maxFlagFraction', maxFlagFraction, ...
              'scaleDays', days <= optimizeIdx(end), 'protected', iborContracts(K, days));
fprintf('Model %s, %d days from %s, %d contracts, threshold %g\n', cmpPeriods.modelFile, ...
        numel(days), datestr(floor(times(days(1))), 'yyyy-mm-dd'), K.nu, threshold);

runs = struct('name', {}, 'flag', {}, 'r', {});
fprintf('\nA: start parameters, robust scale...\n');
[flag, r] = flagQuotes(K, xStart, days, setfield(opts, 'scale', 'robust')); %#ok<SFLD>
runs(end+1) = struct('name', 'A start/robust', 'flag', flag, 'r', r);

resultPath = dataPath('kalmanOptimizeResult.mat');
if isfile(resultPath)
  R = load(resultPath);
  if strcmp(R.modelFile, cmpPeriods.modelFile) && numel(R.xOpt) == numel(xStart)
    fprintf('B: optimized parameters (saved %s), robust scale...\n', R.savedAt);
    [flag, r] = flagQuotes(K, R.xOpt, days, setfield(opts, 'scale', 'robust')); %#ok<SFLD>
    runs(end+1) = struct('name', 'B opt/robust', 'flag', flag, 'r', r);
    fprintf('C: optimized parameters, absolute scale...\n');
    [flag, r] = flagQuotes(K, R.xOpt, days, setfield(opts, 'scale', 'absolute')); %#ok<SFLD>
    runs(end+1) = struct('name', 'C opt/absolute', 'flag', flag, 'r', r);
  else
    fprintf('kalmanOptimizeResult.mat is for %s, not this model: only A is shown.\n', R.modelFile);
  end
else
  fprintf('No kalmanOptimizeResult.mat: only A is shown.\n');
end

%% Summary per run
nObs = sum(runs(1).r.nz);
fprintf('\n%-16s %8s %9s %7s %12s %12s\n', 'Run', 'Flagged', 'Share %', 'Capped', 'Spread med', 'Spread max');
for k = 1:numel(runs)
  rr = runs(k).r;
  rawStd = 1.4826*median(abs(rr.rLOO(~isnan(rr.rLOO)) - median(rr.rLOO(~isnan(rr.rLOO)))));
  fprintf('%-16s %8d %9.3f %7d %12.3g %12.3g\n', runs(k).name, nnz(runs(k).flag), ...
          100*nnz(runs(k).flag)/nObs, nnz(rr.capped), median(rr.rScale), max(rr.rScale));
  fprintf('%-16s robust std of rLOO with these parameters: %.3g (about 1 = calibrated)\n', '', rawStd);
end

%% Overlap with A
if numel(runs) > 1
  fprintf('\nOverlap with A (Jaccard = both / either):\n');
  fprintf('%-16s %8s %8s %8s %8s\n', 'Run', 'Both', 'Only A', 'Only it', 'Jaccard');
  a = runs(1).flag;
  for k = 2:numel(runs)
    b = runs(k).flag;
    both = nnz(a & b); onlyA = nnz(a & ~b); onlyB = nnz(~a & b);
    fprintf('%-16s %8d %8d %8d %8.2f\n', runs(k).name, both, onlyA, onlyB, both/max(1, both+onlyA+onlyB));
  end

  % Per contract, A against the best reference (B)
  b = runs(2).flag;
  fprintf('\nPer contract (A / %s / both):\n', runs(2).name);
  for i = find(any(a | b, 1))
    fprintf('  %-14s %5d %5d %5d\n', ricAll{i}, nnz(a(:,i)), nnz(b(:,i)), nnz(a(:,i) & b(:,i)));
  end

  %% Plots of disagreements
  [ti, ci] = find(xor(a, b));
  if ~isempty(ti) && nExamplePlots > 0
    pick = round(linspace(1, numel(ti), min(nExamplePlots, numel(ti))));
    figure('Name', 'compareFlags: A and B disagree');
    for p = 1:numel(pick)
      t = ti(pick(p)); i = ci(pick(p));
      w = max(days(1), t-15):min(days(end), t+15);
      rows = w - days(1) + 1;
      subplot(ceil(numel(pick)/2), 2, p);
      plot(floor(times(w)), runs(1).r.z(rows, i), '.-', floor(times(w)), runs(2).r.pred(rows, i), '-');
      hold on;
      plot(floor(times(t)), runs(1).r.z(t-days(1)+1, i), 'o', 'MarkerSize', 10);
      hold off;
      datetick('x', 'yyyy-mm-dd', 'keeplimits');
      title(sprintf('%s %s: A %d, B %d', ricAll{i}, datestr(floor(times(t)), 'yyyy-mm-dd'), a(t,i), b(t,i)), ...
            'Interpreter', 'none');
    end
    legend('quote', 'prediction (B)', 'disagreement', 'Location', 'best');
  end
end

function protected = iborContracts(K, days)
% The overnight fixing is never removed (as in kalmanOptimizeMain).
  protected = false(1, K.nu);
  for t = days
    [row, col] = find(K.Iz{t});
    protected(col(K.atAll{t}(row) == K.pl.atIBOR)) = true;
  end
end
