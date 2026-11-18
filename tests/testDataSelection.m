function tests = testDataSelection
% Tests for selectDays and subsetKalmanModel.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

%% selectDays

function testDateWindowIgnoresTimeOfDay(testCase)
times = weekdaysFrom(datenum(2018,1,1), 60) + 22/24;
idx = selectDays(times, datenum(2018,1,3), datenum(2018,1,10));
verifyEqual(testCase, floor(times(idx))', datenum(2018,1,[3 4 5 8 9 10]));
end

function testEmptyBoundsKeepEverything(testCase)
times = weekdaysFrom(datenum(2018,1,1), 30);
verifyEqual(testCase, selectDays(times, [], []), 1:30);
verifyEqual(testCase, selectDays(times, datenum(2018,1,3), []), 3:30);
end

function testNoDaysInWindow(testCase)
times = weekdaysFrom(datenum(2018,1,1), 30);
verifyError(testCase, @() selectDays(times, datenum(2019,1,1), []), 'kalmanrates:NoObservations');
verifyError(testCase, @() selectDays(times, datenum(2018,2,1), datenum(2018,1,1)), 'kalmanrates:NoObservations');
end

%% subsetKalmanModel

function testDaysAreCopiedAndOriginalUnchanged(testCase)
[K, times] = modelFixture();
before = snapshot(K);
idx = selectDays(times, datenum(2018,1,10), datenum(2018,2,20));
KS = subsetKalmanModel(K, idx);
verifyEqual(testCase, snapshot(K), before);
verifyTrue(testCase, isstruct(KS));
verifyEqual(testCase, KS.nt, numel(idx));
verifyEqual(testCase, KS.z, K.z(idx));
verifyEqual(testCase, KS.B, K.B(idx));
verifyEqual(testCase, KS.oTAll, K.oTAll(idx));
verifyEqual(testCase, KS.nz, K.nz(idx));
verifyEqual(testCase, KS.thetaF, K.thetaF);
verifyEqual(testCase, KS.Sigmaw, K.Sigmaw);
end

function testAllDaysGivesSameModel(testCase)
[K, ~] = modelFixture();
KS = subsetKalmanModel(K, 1:K.nt);
before = snapshot(K);
names = fieldnames(before);
for i = 1:numel(names)
  verifyEqual(testCase, KS.(names{i}), before.(names{i}), ['Field ' names{i}]);
end
end

function testInvalidIndices(testCase)
[K, ~] = modelFixture();
verifyError(testCase, @() subsetKalmanModel(K, [3 2]), 'kalmanrates:InvalidIndices');
verifyError(testCase, @() subsetKalmanModel(K, [1 3]), 'kalmanrates:InvalidIndices');
verifyError(testCase, @() subsetKalmanModel(K, [K.nt K.nt+1]), 'kalmanrates:InvalidIndices');
verifyError(testCase, @() subsetKalmanModel(K, []), 'kalmanrates:InvalidIndices');
end

%% Helpers

function times = weekdaysFrom(first, n)
days = first + (0:2*n)';
days = days(weekday(days) >= 2 & weekday(days) <= 6);
times = days(1:n);
end

function [K, times] = modelFixture()
nt = 60;
np = 6;
nc = 16;
nu = 3;
nx = np+nc+nu;
times = weekdaysFrom(datenum(2018,1,1), nt) + 22/24;
indInstrAll = repmat({1:nu}, nt, 1);
prices = arrayfun(@(t) [t; t+0.1; t+0.2], (1:nt)', 'UniformOutput', false);
atAll = repmat({[20; 24; 24]}, nt, 1);
oInd = repmat({[1; 2; 4; 6]}, nt, 1);
oT = arrayfun(@(t) t*10 + (1:5)', (1:nt)', 'UniformOutput', false);
tc = arrayfun(@(t) t + (1:5)'/10, (1:nt)', 'UniformOutput', false);
o = arrayfun(@(t) t + (1:5)'*[1 2], (1:nt)', 'UniformOutput', false);
dummy = repmat({zeros(nc,1)}, nt, 1);
stepDates = [datenum(2018,1,18); datenum(2018,1,25); datenum(2018,3,1) + (0:nc)'*42];
K = clsKalmanIr(times, zeros(nt,np), indInstrAll, true(1,nu), prices, ...
    struct(), atAll, o, oT, oInd, tc, nc, dummy, stepDates, ...
    (1:nx)'/nx, ones(nx,1), (1:nu)', zeros(nx,1), ...
    zeros(nx,1), ones(nx,1), eye(25,np));
end

function saved = snapshot(K)
names = properties(K);
saved = struct();
for i = 1:numel(names)
  saved.(names{i}) = K.(names{i});
end
end
