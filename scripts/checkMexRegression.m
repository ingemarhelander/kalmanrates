function report = checkMexRegression(mode, referenceFile)
% Save old-MEX results before rebuilding, then verify the new MEX against them.
% Run on the full-data computer: this explicitly checks 100/250/500/all dates.
%   checkMexRegression('baseline')
%   compileMex
%   checkMexRegression('verify')
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
if nargin < 1
  mode = 'verify';
end
if nargin < 2
  referenceFile = dataPath('adCountReference.mat');
end
sourceFile = dataPath('kalmanDebug.mat');
root = fileparts(fileparts(mfilename('fullpath')));
expectedMex = fullfile(root, 'mex', ['mexKalmanLogLikelihood.' mexext]);
assert(strcmpi(which('mexKalmanLogLikelihood'), expectedMex), ...
    'kalmanrates:WrongMex', ...
    'MATLAB must use %s. Check which mexKalmanLogLikelihood -all.', expectedMex);
sourceInfo = dir(sourceFile);
binaryInfo = dir(expectedMex);
assert(~isempty(sourceInfo) && ~isempty(binaryInfo), ...
    'kalmanrates:MissingRegressionInput', 'The saved model and MEX binary are required.');
sourceStamp = [sourceInfo.bytes sourceInfo.datenum];
binaryStamp = [binaryInfo.bytes binaryInfo.datenum];

switch mode
  case 'baseline'
    assert(exist(referenceFile, 'file') == 0, ...
        'kalmanrates:ReferenceExists', ...
        'Reference already exists. Use a different referenceFile to preserve it.');
    saved = load(sourceFile, 'K');
    optima = load(dataPath('kalmanXoptAll.mat'), 'xOptAll');
    parameters = optima.xOptAll(:,1);
    assert(saved.K.nt >= 500, 'kalmanrates:InsufficientRegressionHistory', ...
        'At least 500 dates are required to distinguish the sample sizes.');
    clear saved optima
    reference.version = 1;
    reference.parameters = parameters;
    reference.sourceStamp = sourceStamp;
    reference.binaryStamp = binaryStamp;
    reference.limits = [100 250 500 Inf];
    reference.results = cell(size(reference.limits));
    for i=1:numel(reference.limits)
      reference.results{i} = evaluateCase(sourceFile, parameters, reference.limits(i), false);
    end
    save(referenceFile, 'reference');
    report = reference;
    fprintf('Reference saved to %s. Rebuild with compileMex, then run verify.\n', referenceFile);

  case 'verify'
    saved = load(referenceFile, 'reference');
    reference = saved.reference;
    assert(reference.version == 1, 'kalmanrates:ReferenceVersion', ...
        'Unsupported regression reference version.');
    assert(isequal(sourceStamp, reference.sourceStamp), ...
        'kalmanrates:RegressionInputChanged', ...
        'kalmanDebug.mat changed since the baseline was recorded.');
    assert(~isequal(binaryStamp, reference.binaryStamp), ...
        'kalmanrates:MexNotRebuilt', 'Rebuild the MEX binary before verification.');
    testMexDir = fullfile(root, 'tests', 'mex');
    originalPath = path;
    pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
    addpath(testMexDir);
    assert(exist('testAdCountsMex', 'file') == 3, ...
        'kalmanrates:MissingAdTestMex', 'Run compileMex to build the AD test helper.');
    assert(testAdCountsMex() == 6);
    report = cell(size(reference.limits));
    for i=1:numel(reference.limits)
      current = evaluateCase(sourceFile, reference.parameters, reference.limits(i), true);
      previous = reference.results{i};
      assert(isequal(current.times, previous.times) && current.nt == previous.nt, ...
          'kalmanrates:RegressionDateMismatch', 'Observation dates differ from the baseline.');
      assertClose(current.f, previous.f, 'baseline objective');
      assertClose(current.gradient, previous.gradient, 'baseline gradient');
      report{i} = current;
    end
    fprintf('All sample sizes match the baseline; the 250-date gradient checks passed.\n');

  otherwise
    error('kalmanrates:RegressionMode', 'Use baseline or verify.');
end
end

function result = evaluateCase(sourceFile, parameters, limit, verifyNewMex)
saved = load(sourceFile, 'K', 'KI');
[K, ~] = limitKalmanObservations(saved.K, limit, saved.KI.times);
result.nt = K.nt;
result.times = saved.KI.times(1:K.nt);
clear saved
blockSizes = [numel(K.thetaF), numel(K.Sigmav), numel(K.Sigmaw), ...
              numel(K.ax), numel(K.Sigma0), numel(K.Ef)];
assert(numel(parameters) == sum(blockSizes), ...
    'kalmanrates:ParameterDimensions', 'Saved optimum does not match the saved model.');
edges = [0 cumsum(blockSizes)];
K.thetaF = parameters(edges(1)+1:edges(2));
K.Sigmav = parameters(edges(2)+1:edges(3));
K.Sigmaw = parameters(edges(3)+1:edges(4));
K.ax = parameters(edges(4)+1:edges(5));
K.Sigma0 = parameters(edges(5)+1:edges(6));
K.Ef(:) = parameters(edges(6)+1:edges(7));

% Capture the native diagnostic where MATLAB redirects it to the command window.
output = evalc('[f, gradient] = kalmanObj(parameters, K, 0);');
if ~isempty(output)
  fprintf('%s', output);
end
assert(isfinite(f) && all(isfinite(gradient)), ...
    'kalmanrates:NonfiniteRegression', 'Objective and gradient must be finite.');
if verifyNewMex
  assert(isempty(strfind(output, 'AD: Incorrect number')), ...
      'kalmanrates:AdCountMismatch', 'The rebuilt MEX still reports an AD count mismatch.');
end
fOnly = kalmanObj(parameters, K, 0);
assertClose(fOnly, f, 'objective with and without gradient');
result.f = f;
result.gradient = gradient;
if verifyNewMex && limit == 250
  % Exercise both sparse covariance products: Sigmaw in each state block,
  % and Sigmav. Use positive entries so central differences stay feasible.
  groups = {edges(3)+(1:K.np), ...
            edges(3)+(K.np+1:K.ns), ...
            edges(3)+(K.ns+1:numel(K.Sigmaw)), ...
            edges(2)+(1:numel(K.Sigmav))};
  for j=1:numel(groups)
    group = groups{j};
    positive = find(parameters(group) > 0, 1);
    assert(~isempty(positive), 'kalmanrates:FiniteDifferenceParameters', ...
        'Each tested covariance block needs at least one positive parameter.');
    index = group(positive);
    checkCentralDifference(parameters, K, index, gradient(index));
  end
end
fprintf('Checked %d observations (limit %g).\n', result.nt, limit);
end

function checkCentralDifference(parameters, K, index, analytic)
relativeSteps = [1e-3 1e-4 1e-5];
bestError = Inf;
for relativeStep = relativeSteps
  h = min(relativeStep*max(abs(parameters(index)), 1e-6), parameters(index)/4);
  plus = parameters;
  minus = parameters;
  plus(index) = plus(index)+h;
  minus(index) = minus(index)-h;
  fPlus = kalmanObj(plus, K, 0);
  fMinus = kalmanObj(minus, K, 0);
  numeric = (fPlus-fMinus)/(2*h);
  if isfinite(numeric)
    scaledError = abs(numeric-analytic)/(1e-5+1e-3*max(abs(numeric),abs(analytic)));
    bestError = min(bestError, scaledError);
  end
end
assert(bestError <= 1, 'kalmanrates:GradientMismatch', ...
    'Central differences disagree with gradient component %d (scaled error %g).', index, bestError);
end

function assertClose(actual, expected, label)
assert(isequal(size(actual), size(expected)) && ...
       all(isfinite(actual(:))) && all(isfinite(expected(:))) && ...
       all(abs(actual(:)-expected(:)) <= 1e-8+1e-10*abs(expected(:))), ...
    'kalmanrates:RegressionMismatch', 'Mismatch in %s.', label);
end
