function tests = testDataLimit
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

function testPrefixAndOutputNames(testCase)
times = datenum(2005,8,15) + (0:499)';
[indices, suffix] = selectObservations(numel(times), 250, times);
verifyEqual(testCase, indices, 1:250);
verifyEqual(testCase, suffix, '_sample250');
for limit = [500, 1000, Inf]
  [indices, suffix] = selectObservations(numel(times), limit, times);
  verifyEqual(testCase, indices, 1:500);
  verifyEqual(testCase, suffix, '');
end
end

function testInvalidLimitsAndDates(testCase)
invalid = {0, -1, 1.5, NaN, -Inf, [], [100 250], '250', 1+1i, true};
for i=1:numel(invalid)
  verifyError(testCase, @() selectObservations(500, invalid{i}), ...
              'kalmanrates:InvalidObservationLimit');
end
verifyError(testCase, @() selectObservations(0, 250), 'kalmanrates:NoObservations');
verifyError(testCase, @() selectObservations(500, 250, 1:499), ...
            'kalmanrates:DateAlignment');
end

function testSavedModelAlignmentAndParameters(testCase)
[K, times] = modelFixture();
original = snapshot(K);
[K, suffix] = limitKalmanObservations(K, 250, times);
verifyEqual(testCase, K.nt, 250);
verifyEqual(testCase, suffix, '_sample250');
names = fieldnames(original);
for i=1:numel(names)
  name = names{i};
  expected = original.(name);
  if iscell(expected) || strcmp(name, 'nz')
    expected = expected(1:250);
  elseif strcmp(name, 'nt')
    expected = 250;
  end
  verifyEqual(testCase, K.(name), expected, ['Incorrect field: ' name]);
end
% The prefix retains the original transition and the first observation prior.
verifyEqual(testCase, K.B{1}, original.B{1});
verifyEqual(testCase, K.x0, original.x0);
end

function testFullModelUnchanged(testCase)
for limit = [1000, Inf]
  [K, times] = modelFixture();
  original = snapshot(K);
  [K, suffix] = limitKalmanObservations(K, limit, times);
  verifyEqual(testCase, snapshot(K), original);
  verifyEqual(testCase, suffix, '');
end
end

function testMisalignedModelIsNotPartiallyTrimmed(testCase)
[K, times] = modelFixture();
K.oTAll = K.oTAll(1:end-1);
original = snapshot(K);
verifyError(testCase, @() limitKalmanObservations(K, 250, times), ...
            'kalmanrates:ModelAlignment');
verifyEqual(testCase, snapshot(K), original);
end

function testPcaWithoutOutOfSampleDates(testCase)
previousRng = rng;
rngCleanup = onCleanup(@() rng(previousRng)); %#ok<NASGU>
rng(42);
previousVisibility = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
visibilityCleanup = onCleanup(@() set(groot, 'defaultFigureVisible', previousVisibility)); %#ok<NASGU>
oldFigures = findall(groot, 'Type', 'figure');
figureCleanup = onCleanup(@() close(setdiff(findall(groot, 'Type', 'figure'), oldFigures))); %#ok<NASGU>
times = datenum(2005,8,15) + (0:249)';
firstDates = floor(times);
lastDates = firstDates + 20;
fH = cumsum(randn(250,25)*1e-4,1);
piH = cumsum(randn(250,25)*1e-5,1);
[fH, piH, lastDates, nIn, ef, Ef] = pcaIRS(times, fH, piH, ...
    firstDates, lastDates, datenum(2020,1,1));
verifyEqual(testCase, nIn, 250);
verifyEqual(testCase, size(Ef), [25 6]);
verifyTrue(testCase, all(isfinite([fH(:); piH(:); ef(:); Ef(:)])));
verifyEqual(testCase, lastDates-firstDates, ones(250,1)*25);
verifyError(testCase, @() pcaIRS(times(1:7), fH(1:7,:), piH(1:7,:), ...
    firstDates(1:7), lastDates(1:7), datenum(2020,1,1)), ...
    'kalmanrates:InsufficientPcaHistory');
end

function [K, times] = modelFixture()
nt = 500;
np = 6;
nc = 16;
nu = 3;
nx = np+nc+nu;
times = datenum(2005,8,15) + (0:nt-1)';
indInstrAll = repmat({1:nu},nt,1);
prices = arrayfun(@(t) ones(nu,1)*t, (1:nt)', 'UniformOutput', false);
dummy = repmat({zeros(nu,1)},nt,1);
stepDates = times(end)+(1:nc)'*30;
K = clsKalmanIr(times, zeros(nt,np), indInstrAll, true(1,nu), prices, ...
    struct(), dummy, dummy, dummy, dummy, dummy, nc, dummy, stepDates, ...
    ones(nx,1), ones(nx,1), ones(nu,1), zeros(nx,1), ...
    (1:nx)', ones(nx,1), eye(25,np));
end

function saved = snapshot(K)
names = properties(K);
saved = struct();
for i=1:numel(names)
  saved.(names{i}) = K.(names{i});
end
end
