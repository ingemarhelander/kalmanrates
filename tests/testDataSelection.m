function tests = testDataSelection
% Tests for selectDays, selectInstruments and subsetKalmanModel.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

%% selectDays

function testDateWindowIgnoresTimeOfDay(testCase)
times = weekdaysFrom(datenum(2018,1,1), 60) + 22/24;
idx = selectDays(times, struct('startDate', '2018-01-03', 'endDate', '2018-01-10'));
verifyEqual(testCase, floor(times(idx))', datenum(2018,1,[3 4 5 8 9 10]));
idx = selectDays(times, struct('startDate', datenum(2018,1,3), 'endDate', datenum(2018,1,10)));
verifyEqual(testCase, floor(times(idx))', datenum(2018,1,[3 4 5 8 9 10]));
end

function testDefaultsKeepEverything(testCase)
times = weekdaysFrom(datenum(2018,1,1), 30);
verifyEqual(testCase, selectDays(times, struct()), 1:30);
verifyEqual(testCase, selectDays(times), 1:30);
end

function testWeekdayFilter(testCase)
times = weekdaysFrom(datenum(2018,1,1), 60);
idx = selectDays(times, struct('weekdays', 'wed'));
verifyTrue(testCase, all(weekday(times(idx)) == 4));
verifyEqual(testCase, numel(idx), 12);
idx = selectDays(times, struct('weekdays', 'Mon, thursday'));
verifyEqual(testCase, unique(weekday(times(idx)))', [2 5]);
idx = selectDays(times, struct('weekdays', {{'fri'}}));
verifyTrue(testCase, all(weekday(times(idx)) == 6));
end

function testDayStepAfterOtherFilters(testCase)
times = weekdaysFrom(datenum(2018,1,1), 60);
idx = selectDays(times, struct('dayStep', 5));
verifyEqual(testCase, idx, 1:5:60);
idx = selectDays(times, struct('startDate', '2018-01-02', 'weekdays', 'mon,wed', 'dayStep', 2));
allIdx = selectDays(times, struct('startDate', '2018-01-02', 'weekdays', 'mon,wed'));
verifyEqual(testCase, idx, allIdx(1:2:end));
end

function testInvalidDaySelections(testCase)
times = weekdaysFrom(datenum(2018,1,1), 30);
verifyError(testCase, @() selectDays(times, struct('startDate', '2018/01/05')), 'kalmanrates:InvalidDate');
verifyError(testCase, @() selectDays(times, struct('weekdays', 'xyz')), 'kalmanrates:InvalidWeekday');
verifyError(testCase, @() selectDays(times, struct('weekdays', 'sat')), 'kalmanrates:NoObservations');
verifyError(testCase, @() selectDays(times, struct('startDate', '2019-01-01')), 'kalmanrates:NoObservations');
verifyError(testCase, @() selectDays(times, struct('startDate', '2018-02-01', 'endDate', '2018-01-01')), ...
            'kalmanrates:InvalidDateRange');
for bad = {0, 1.5, -2, NaN, [1 2]}
  verifyError(testCase, @() selectDays(times, struct('dayStep', bad{1})), 'kalmanrates:InvalidDayStep');
end
end

%% selectInstruments

function testTenorNames(testCase)
ric = {'EONIA=', 'EUREONSW=', 'EUREON3M=', 'EUREON10Y=', 'EUREON11Y='};
[keep, labels] = selectInstruments(ric, 4, 'all');
verifyEqual(testCase, keep, true(1,4));
verifyEqual(testCase, labels, {'EONIA', '1W', '3M', '10Y'});
verifyEqual(testCase, selectInstruments(ric, 4, '1w, 10Y'), logical([0 1 0 1]));
verifyEqual(testCase, selectInstruments(ric, 4, {'EUREON3M=', 'eonia'}), logical([1 0 1 0]));
verifyEqual(testCase, selectInstruments(ric, 4, 'SW'), logical([0 1 0 0]));
verifyError(testCase, @() selectInstruments(ric, 4, '11Y'), 'kalmanrates:UnknownTenor');
verifyError(testCase, @() selectInstruments(ric(1:3), 4, 'all'), 'kalmanrates:ModelAlignment');
end

%% subsetKalmanModel

function testDaysAreCopiedAndOriginalUnchanged(testCase)
[K, times] = modelFixture();
before = snapshot(K);
idx = selectDays(times, struct('startDate', '2018-01-10', 'endDate', '2018-02-20'));
KS = subsetKalmanModel(K, idx);
verifyEqual(testCase, snapshot(K), before);
verifyTrue(testCase, isstruct(KS));
verifyEqual(testCase, KS.nt, numel(idx));
verifyEqual(testCase, KS.z, K.z(idx));
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

function testCentralBankStepMovesToNextKeptDay(testCase)
[K, times] = modelFixture();
shift = find(isShiftDay(K));
verifyNotEmpty(testCase, shift);
idx = selectDays(times, struct('weekdays', 'wed'));
KS = subsetKalmanModel(K, idx);
% Every original step lands on exactly one kept day: the first kept day at or after it.
expected = false(1, numel(idx));
for s = shift(shift > idx(1))
  expected(find(idx >= s, 1)) = true;
end
verifyEqual(testCase, isShiftDay(KS), expected);
verifyEqual(testCase, sum(isShiftDay(KS)), sum(shift > idx(1)));
for j = find(expected)
  verifyEqual(testCase, KS.B{j}, K.B{shift(1)});
  verifyEqual(testCase, KS.F{j}, sparse(KS.A{j}*diag(KS.thetaF)*KS.B{j}));
end
end

function testTooManyStepsInOneGap(testCase)
[K, ~] = modelFixture();
shift = find(isShiftDay(K));
verifyError(testCase, @() subsetKalmanModel(K, [1, shift(2)+1]), 'kalmanrates:TooManySteps');
end

function testDropContracts(testCase)
[K, ~] = modelFixture();
% Instruments: 1 = IBOR (one row in o), 2 and 3 = OIS (two rows each).
KS = subsetKalmanModel(K, 1:K.nt, logical([1 0 1]));
t = 7;
verifyEqual(testCase, KS.nz(t), 2);
verifyEqual(testCase, KS.z{t}, K.z{t}([1 3]));
verifyEqual(testCase, KS.atAll{t}, K.atAll{t}([1 3]));
verifyEqual(testCase, KS.oTAll{t}, K.oTAll{t}([1 4 5]));
verifyEqual(testCase, KS.tcAll{t}, K.tcAll{t}([1 4 5]));
verifyEqual(testCase, KS.oAll{t}, K.oAll{t}([1 4 5], :));
verifyEqual(testCase, KS.oIndAll{t}, [1; 2; 4]);
verifyEqual(testCase, full(KS.Iz{t}), [1 0 0; 0 0 1]);
verifyEqual(testCase, full(KS.G{t}), [1 0 0; 0 0 1]);
verifyEqual(testCase, KS.Cv{t}, K.Cv{t}([1 3]));
verifyEqual(testCase, KS.nuv{t}, K.nuv{t}([1 3]));
% Parameter and state sizes do not change.
verifyEqual(testCase, KS.nu, K.nu);
verifyEqual(testCase, numel(KS.Sigmav), numel(K.Sigmav));
verifyEqual(testCase, KS.nx, K.nx);
% Dropping only the OIS contracts keeps one row per day.
KS = subsetKalmanModel(K, 1:K.nt, logical([1 0 0]));
verifyEqual(testCase, KS.nz, ones(K.nt, 1));
verifyError(testCase, @() subsetKalmanModel(K, 1:K.nt, false(1,3)), 'kalmanrates:NoInstruments');
end

function testInvalidIndices(testCase)
[K, ~] = modelFixture();
verifyError(testCase, @() subsetKalmanModel(K, [3 2]), 'kalmanrates:InvalidIndices');
verifyError(testCase, @() subsetKalmanModel(K, [1 K.nt+1]), 'kalmanrates:InvalidIndices');
verifyError(testCase, @() subsetKalmanModel(K, []), 'kalmanrates:InvalidIndices');
end

%% Helpers

function times = weekdaysFrom(first, n)
days = first + (0:2*n)';
days = days(weekday(days) >= 2 & weekday(days) <= 6);
times = days(1:n);
end

function shift = isShiftDay(K)
shift = false(1, K.nt);
for t = 1:K.nt
  shift(t) = nnz(K.B{t} - speye(K.nx)) > 0;
end
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
% Central bank steps on two Thursdays, one week apart, and later dates.
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
