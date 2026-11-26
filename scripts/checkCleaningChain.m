% checkCleaningChain.m
% End-to-end check of the cleaning on a COPY of the project, so data/ in the
% project is never touched:
%   cleanMarketData -> buildKalmanModel -> comparePasses -> useDataSet
%
% The copy is made in tempdir (kalmanrates_cleaning_check). It needs
% data/kalmanOptimizeResult.mat, the model file it refers to, EUR.mat and
% EUR100.mat (or EUR_raw.mat and EUR100_raw.mat) and the MEX files.
%
% The pass-1 parameters may not be calibrated (then nothing can be flagged),
% so the copy uses the same parameters with the noise rescaled and the
% calibration check turned off. The flags only test the mechanics; they are
% not a real cleaning.
%
% Takes a few minutes: EUR100.mat is copied and the model is built once.

project = fileparts(fileparts(mfilename('fullpath')));
SB = fullfile(tempdir, 'kalmanrates_cleaning_check');

%% 0. Fresh copy of the project
fprintf('Copying the project to %s ...\n', SB);
if isfolder(SB)
  warning('off', 'MATLAB:rmpath:DirNotFound');
  rmpath(genpath(SB));
  rmdir(SB, 's');
end
mkdir(fullfile(SB, 'data', 'raw'));
for d = {'measurement', 'model', 'scripts', 'mex', 'tests'}
  copyfile(fullfile(project, d{1}), fullfile(SB, d{1}));
end
copyfile(fullfile(project, '*.m'), SB);
copyfile(fullfile(project, 'data', 'raw', 'datesEURECB.xlsx'), fullfile(SB, 'data', 'raw'));
source = fullfile(project, 'data');
Rreal = load(fullfile(source, 'kalmanOptimizeResult.mat'));
if isfile(fullfile(source, 'EUR_raw.mat'))  % the raw data, whichever set is active
  copyfile(fullfile(source, 'EUR_raw.mat'), fullfile(SB, 'data', 'EUR.mat'));
  copyfile(fullfile(source, 'EUR100_raw.mat'), fullfile(SB, 'data', 'EUR100.mat'));
else
  copyfile(fullfile(source, 'EUR.mat'), fullfile(SB, 'data'));
  copyfile(fullfile(source, 'EUR100.mat'), fullfile(SB, 'data'));
end
copyfile(fullfile(source, 'kalmanOptimizeResult.mat'), fullfile(SB, 'data'));
copyfile(fullfile(source, Rreal.modelFile), fullfile(SB, 'data'));
% Allow writing although the calibration check fails (only in the copy)
f = fullfile(SB, 'scripts', 'cleanMarketData.m');
txt = fileread(f);
assert(contains(txt, 'allowUncalibrated = false;'), 'allowUncalibrated not found in cleanMarketData.m');
txt = strrep(txt, 'allowUncalibrated = false;', 'allowUncalibrated = true; ');
fid = fopen(f, 'w'); fwrite(fid, txt); fclose(fid);

% Work in the copy: the current folder comes before the path, so cd is needed
% for dataPath to point to the copy. Folder and path are restored at the end,
% also when a check fails.
clear functions  % forget the functions of the project
cd(SB);
run(fullfile(SB, 'setupPaths.m'));
try
assert(startsWith(dataPath('x'), SB), 'dataPath does not point to the copy');

%% 1. Pass-1 result with rescaled noise
R = load(dataPath('kalmanOptimizeResult.mat'));
[K0, times0, ric0] = loadKalmanModel(R.modelFile);
n1 = numel(K0.thetaF); n2 = numel(K0.Sigmav); n3 = numel(K0.Sigmaw);
iv = n1+(1:n2); iw = n1+n2+(1:n3);
R.xOpt(iv) = 1e-12;
R.xOpt(iw(1:K0.np)) = R.xOpt(iw(1:K0.np))*0.000255;
R.xOpt(iw(K0.np+1:K0.ns)) = R.xOpt(iw(K0.np+1:K0.ns))*9.56e-07;
R.xOpt(iw(K0.ns+1:end)) = R.xOpt(iw(K0.ns+1:end))*2.71e-08;
R.savedAt = 'check pass 1';
save(dataPath('kalmanOptimizeResult.mat'), '-struct', 'R');
nz0 = sum(K0.nz);

%% 2. Clean
clear dryRun
run(fullfile(SB, 'scripts', 'cleanMarketData.m'));
close all
for f = {'EUR_raw.mat', 'EUR_clean.mat', 'EUR100_raw.mat', 'dataSet.mat', 'noisyQuotes.mat', 'kalmanOptimizeResult_pass1.mat'}
  assert(isfile(dataPath(f{1})), 'missing %s', f{1});
end
Q = load(dataPath('noisyQuotes.mat'));
assert(height(Q.flags) > 0, 'nothing was flagged in the check');
cleanData = load(dataPath('EUR_clean.mat')); rawData = load(dataPath('EUR_raw.mat'));
assert(isequaln(load(dataPath('EUR.mat')), cleanData), 'EUR.mat is not the cleaned data');
nNaNnew = nnz(isnan(cleanData.data{1})) - nnz(isnan(rawData.data{1}));
assert(nNaNnew == height(Q.flags), 'BID NaN (%d) does not match the flags (%d)', nNaNnew, height(Q.flags));
assert(~any(strcmp(Q.flags.ric, 'EONIA=')), 'protected contract flagged');
assert(strcmp(useDataSet(), 'clean'));
fprintf('\nCHECK 1 OK: %d flags -> %d BID set to NaN, EONIA= not touched\n', height(Q.flags), nNaNnew);

%% 3. Build the model again on the cleaned data
buildKalmanModel;
[K1, times1, ric1] = loadKalmanModel(R.modelFile);
assert(isequal(floor(times0), floor(times1)), 'days changed');
assert(isequal(ric0, ric1), 'contracts changed');
maxDiff = 0;
for t = 1:K0.nt
  [r0, c0] = find(K0.Iz{t}); z0 = nan(1, K0.nu); z0(c0) = K0.z{t}(r0);
  [r1, c1] = find(K1.Iz{t}); z1 = nan(1, K1.nu); z1(c1) = K1.z{t}(r1);
  both = ~isnan(z0) & ~isnan(z1);
  maxDiff = max([maxDiff abs(z0(both) - z1(both))]);
  for c = find(~isnan(z0) & isnan(z1))
    assert(any(Q.flags.date == floor(times0(t)) & strcmp(Q.flags.ric, ric0{c})), ...
           'quote %s %s disappeared without a flag', datestr(times0(t)), ric0{c});
  end
end
assert(nz0 - sum(K1.nz) == height(Q.flags), 'the model lost %d quotes, flags %d', nz0 - sum(K1.nz), height(Q.flags));
assert(maxDiff == 0, 'kept quotes changed');
fprintf('CHECK 2 OK: the model lost exactly the %d flagged quotes, the other %d are unchanged\n', ...
        height(Q.flags), sum(K1.nz));

%% 4. Fake pass 2 and compare
R2 = R; R2.xOpt(iv) = R2.xOpt(iv)*0.5; R2.savedAt = 'check pass 2';
save(dataPath('kalmanOptimizeResult.mat'), '-struct', 'R2');
run(fullfile(SB, 'scripts', 'comparePasses.m'));
fprintf('CHECK 3 OK: comparePasses ran\n');

%% 5. Switching between the data sets
useDataSet('raw');
assert(isequaln(load(dataPath('EUR.mat')), rawData), 'raw not restored');
useDataSet('clean');
assert(isequaln(load(dataPath('EUR.mat')), cleanData), 'clean not restored');
fprintf('CHECK 4 OK: useDataSet switches between raw and clean\n');
checkError = [];
catch checkError
end

%% Back to the project
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(genpath(SB));
clear functions
cd(project);
run(fullfile(project, 'setupPaths.m'));
if ~isempty(checkError)
  fprintf('\nCHECK FAILED (the data in %s was not touched)\n', fullfile(project, 'data'));
  rethrow(checkError);
end
fprintf('\nALL CHECKS OK (the data in %s was not touched)\n', fullfile(project, 'data'));
