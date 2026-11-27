function tests = testPrefilter
% Tests for flagQuotes, removeQuotes and the rScale/burnIn options of
% innovationsLOO, on the synthetic model of testInnovationsLOO (one level
% factor, eight IBOR-like contracts).
% Kör från projektroten med:  setupPaths; runtests('tests/testPrefilter.m')
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

function testRemoveNothingKeepsTheModel(testCase)
[KS, x] = fixture(struct());
KR = removeQuotes(KS, false(KS.nt, KS.nu));
verifyEqual(testCase, KR, KS);
end

function testRemovedModelIsTheFilterWithoutTheQuote(testCase)
% Filtering the model without the flagged quote is the same as the filter
% in innovationsLOO that skips it: the second ("ordinary") run sees exactly
% the data that the first run decided to keep.
[KS, x] = fixture(struct('spike', [20 3 0.003 ; 41 6 -0.004]));
r = innovationsLOO(KS, x);
verifyEqual(testCase, nnz(r.flag), 2);
[KR, ~, n] = removeQuotes(KS, r.flag);
verifyEqual(testCase, n, 2);
verifyEqual(testCase, KR.nz([20 41]), [7 ; 7]);
ref = filterDays(KR, x);
verifyEqual(testCase, ref.f, r.f, 'RelTol', 1e-10);
verifyEqual(testCase, ref.pred(~r.flag), r.pred(~r.flag), 'AbsTol', 1e-14);
verifyTrue(testCase, isnan(ref.z(20, 3)) && isnan(ref.z(41, 6)));
% Nothing is left to flag in the removed model
r2 = innovationsLOO(KR, x);
verifyFalse(testCase, any(r2.flag(:)));
end

function testRobustScaleWorksWithUncalibratedParameters(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003]));
x = scaleNoise(KS, x, 1e3);  % noise variances 1000 times too large
rAbs = innovationsLOO(KS, x);
verifyFalse(testCase, any(rAbs.flag(:)), 'Uncalibrated: the absolute threshold sees nothing');
[flag, r] = flagQuotes(KS, x, 1:KS.nt, struct('burnIn', 5));
verifyEqual(testCase, find(flag), sub2ind(size(flag), 20, 3));
verifyLessThan(testCase, max(r.rScale), 0.2, 'The spread reflects the too large noise');
end

function testRobustScaleGivesNoFlagsOnCleanData(testCase)
[KS, x] = fixture(struct());
flag = flagQuotes(KS, scaleNoise(KS, x, 1e3), 1:KS.nt, struct('burnIn', 5));
verifyFalse(testCase, any(flag(:)));
end

function testAbsoluteScaleIsInnovationsLOO(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003]));
[flag, r] = flagQuotes(KS, x, 1:KS.nt, struct('scale', 'absolute', 'burnIn', 0));
ref = innovationsLOO(KS, x);
verifyEqual(testCase, flag, ref.flag);
verifyEqual(testCase, r.rScale, ones(1, KS.nu));
end

function testBurnInDaysAreNotFlagged(testCase)
[KS, x] = fixture(struct('spike', [3 2 0.003]));
r = innovationsLOO(KS, x, struct('burnIn', 5));
verifyFalse(testCase, any(r.flag(:)));
verifyGreaterThan(testCase, abs(r.rLOO(3, 2)), 4, 'rLOO is still computed during burn-in');
r = innovationsLOO(KS, x);
verifyTrue(testCase, r.flag(3, 2));
end

function testFlagsOnSubsetOfDays(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003 ; 50 4 0.003]));
flag = flagQuotes(KS, x, 30:KS.nt, struct('scale', 'absolute', 'burnIn', 0));
verifySize(testCase, flag, [KS.nt KS.nu]);
verifyEqual(testCase, find(flag), sub2ind(size(flag), 50, 4), 'Day 20 is outside the days');
end

function testScaleDaysAreUsedForTheSpread(testCase)
[KS, x] = fixture(struct());
x = scaleNoise(KS, x, 1e3);
scaleDays = (1:KS.nt) <= 40;
[~, r1] = flagQuotes(KS, x, 1:KS.nt, struct('burnIn', 5, 'scaleDays', scaleDays, 'minQuotes', 10));
% Make the last days much noisier: the spread from days 1..40 must not change
KN = KS;
for t = 41:KS.nt
  KN.z{t} = KN.z{t} + 5e-4*randn(KS.nu, 1);
end
[~, r2] = flagQuotes(KN, x, 1:KS.nt, struct('burnIn', 5, 'scaleDays', scaleDays, 'minQuotes', 10));
verifyEqual(testCase, r2.rScale, r1.rScale, 'RelTol', 1e-12);
end

function testPricingRowsFollowTheQuote(testCase)
% OIS quotes have two pricing rows, IBOR one: removing a quote must drop
% its own rows and renumber oIndAll.
S = struct('nt', 1, 'nu', 3);
S.Iz = {speye(3)}; S.G = {speye(3)};
S.z = {[1 ; 2 ; 3]}; S.atAll = {[1 ; 2 ; 2]}; S.nz = 3;
S.oIndAll = {[1 ; 2 ; 4 ; 6]};
S.oTAll = {[10 ; 20 ; 21 ; 30 ; 31]};
S.tcAll = {[0.1 ; NaN ; 0.2 ; NaN ; 0.3]};
S.oAll = {(1:5)' * [1 1]};
KR = removeQuotes(S, [false true false]);
verifyEqual(testCase, KR.z{1}, [1 ; 3]);
verifyEqual(testCase, KR.atAll{1}, [1 ; 2]);
verifyEqual(testCase, full(KR.Iz{1}), [1 0 0 ; 0 0 1]);
verifyEqual(testCase, KR.oIndAll{1}, [1 ; 2 ; 4]);
verifyEqual(testCase, KR.oTAll{1}, [10 ; 30 ; 31]);
verifyEqual(testCase, KR.tcAll{1}, [0.1 ; NaN ; 0.3]);
verifyEqual(testCase, KR.oAll{1}, [1 1 ; 4 4 ; 5 5]);
verifyEqual(testCase, KR.nz, 2);
% KI is reduced the same way
KI = struct('indInstrAll', {{[1 2 3]}}, 'priceAll', {{[1 ; 2 ; 3]}});
[~, KIR] = removeQuotes(S, [false true false], KI);
verifyEqual(testCase, KIR.indInstrAll{1}, [1 3]);
verifyEqual(testCase, KIR.priceAll{1}, [1 ; 3]);
end

function testRemovingEveryQuoteOfADayIsAnError(testCase)
[KS, ~] = fixture(struct());
flag = false(KS.nt, KS.nu); flag(10, :) = true;
verifyError(testCase, @() removeQuotes(KS, flag), 'kalmanrates:InvalidFlags');
verifyError(testCase, @() removeQuotes(KS, false(2, 2)), 'kalmanrates:InvalidFlags');
end

%% Helpers

function x = scaleNoise(KS, x, k)
% Multiplies Sigmav and Sigmaw (the noise variances) by k.
nv = numel(KS.Sigmav); nw = numel(KS.Sigmaw); nth = numel(KS.thetaF);
ind = nth + (1:nv+nw);
x(ind) = x(ind) * k;
end

function [KS, x] = fixture(opts)
% Same model as in testInnovationsLOO.
rng(7);
nt = 60; nu = 8; np = 1; nc = 0; ns = np + nc; nx = ns + nu;
T = round(linspace(30, 365, nu))';
nT = max(T) + 1;

KS = struct();
KS.nt = nt; KS.nu = nu; KS.np = np; KS.nc = nc; KS.ns = ns;
KS.pl = struct('atIBOR', 1, 'atOISG', 2);
KS.thetaF = [1 ; 0.5*ones(nu,1)];
KS.Sigmav = 1e-8*ones(nu,1);
KS.Sigmaw = [1e-6 ; 1e-12*ones(nu,1)];
KS.ax = [0.02 ; zeros(nu,1)];
KS.Sigma0 = [1e-6 ; 1e-12*ones(nu,1)];
KS.Ef = ones(nT, np);
x = [KS.thetaF ; KS.Sigmav ; KS.Sigmaw ; KS.ax ; KS.Sigma0 ; KS.Ef(:)];

level = 0.02 + cumsum([0 ; 0.001*randn(nt-1,1)]);
if isfield(opts, 'jump')
  level(opts.jump(1):end) = level(opts.jump(1):end) + opts.jump(2);
end
tc = (T - 1)/365;
g = @(lv) (exp(lv*(T-1)/365) - 1) ./ tc;

[KS.A, KS.B, KS.D, KS.G, KS.Iz, KS.z, KS.atAll, KS.oTAll, KS.oIndAll, KS.tcAll, KS.datesStepRelAll] = ...
  deal(cell(nt,1));
% Fields of clsKalmanIr that the filter does not use, so that
% subsetKalmanModel (used by flagQuotes) accepts the struct
[KS.x, KS.w, KS.v, KS.F, KS.Cw, KS.Cx, KS.Cv, KS.nuv, KS.oAll] = deal(cell(nt,1));
KS.nz = nu*ones(nt,1);
for t = 1:nt
  KS.A{t} = speye(nx); KS.B{t} = speye(nx); KS.D{t} = speye(nx);
  KS.G{t} = speye(nu);
  KS.Iz{t} = speye(nu);
  KS.atAll{t} = ones(nu,1);
  KS.oTAll{t} = T;
  KS.oIndAll{t} = (1:nu+1)';
  KS.tcAll{t} = tc;
  KS.datesStepRelAll{t} = zeros(0,1);
  KS.z{t} = g(level(t)) + 1e-4*randn(nu,1);
end
if isfield(opts, 'spike')
  for k = 1:size(opts.spike, 1)
    s = opts.spike(k,:);
    KS.z{s(1)}(s(2)) = KS.z{s(1)}(s(2)) + s(3);
  end
end
end
