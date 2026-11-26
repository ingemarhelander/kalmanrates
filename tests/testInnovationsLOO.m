function tests = testInnovationsLOO
% Tests for innovationsLOO on a small synthetic model: one level factor and
% eight IBOR-like contracts that all follow the level.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

function testNoThresholdIsSameFilterAsFilterDays(testCase)
[KS, x] = fixture(struct());
r = innovationsLOO(KS, x, struct('threshold', Inf));
ref = filterDays(KS, x);
verifyEqual(testCase, r.f, ref.f, 'RelTol', 1e-10);
verifyEqual(testCase, r.d2, ref.d2, 'RelTol', 1e-10);
verifyEqual(testCase, r.pred, ref.pred, 'AbsTol', 1e-14);
verifyEqual(testCase, r.z, ref.z);
verifyFalse(testCase, any(r.flag(:)));
end

function testCleanDataGivesNoFlags(testCase)
[KS, x] = fixture(struct());
r = innovationsLOO(KS, x);
verifyFalse(testCase, any(r.flag(:)));
q = r.rLOO(~isnan(r.rLOO));
verifyLessThan(testCase, abs(mean(q)), 0.2);
verifyGreaterThan(testCase, std(q), 0.7);
verifyLessThan(testCase, std(q), 1.3);
end

function testSingleSpikeIsFlagged(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003]));
r = innovationsLOO(KS, x);
verifyTrue(testCase, r.flag(20, 3));
verifyEqual(testCase, nnz(r.flag), 1);
verifyGreaterThan(testCase, abs(r.rRemoved(20, 3)), 4);
% The ordinary innovation does not single it out as clearly
verifyGreaterThan(testCase, abs(r.rLOO(20, 3)), 2*abs(r.rStd(20, 3)));
end

function testCommonMoveIsNotFlagged(testCase)
[KS, x] = fixture(struct('jump', [30 0.008]));
r = innovationsLOO(KS, x);
verifyTrue(testCase, all(abs(r.rStd(30, :)) > 4), 'The move should be large for every contract');
verifyFalse(testCase, any(r.flag(:)));
end

function testSpikeOnMoveDayIsFlagged(testCase)
[KS, x] = fixture(struct('jump', [30 0.008], 'spike', [30 5 0.003]));
r = innovationsLOO(KS, x);
verifyEqual(testCase, find(r.flag(30, :)), 5);
verifyEqual(testCase, nnz(r.flag), 1);
end

function testTwoSpikesSameDayAreBothFlagged(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003 ; 20 6 -0.004]));
r = innovationsLOO(KS, x, struct('maxFlagFraction', 0.3));  % 8 contracts: at most 2
verifyEqual(testCase, find(r.flag(20, :)), [3 6]);
verifyEqual(testCase, r.nFlag(20), 2);
verifyFalse(testCase, r.capped(20));
% With the default (20 % of 8 = 1) the day is capped after one removal
r = innovationsLOO(KS, x);
verifyEqual(testCase, r.nFlag(20), 1);
verifyTrue(testCase, r.capped(20));
end

function testRemovedQuoteDoesNotAffectLaterDays(testCase)
[KSa, x] = fixture(struct());
[KSb, ~] = fixture(struct('spike', [20 3 0.003]));
ra = innovationsLOO(KSa, x);
rb = innovationsLOO(KSb, x);
% Without the bad quote day 20 is updated with 7 good quotes, so the next
% days are close to the clean run
verifyEqual(testCase, rb.pred(21:end, :), ra.pred(21:end, :), 'AbsTol', 2e-4);
ref = filterDays(KSb, x);  % keeps the spike
verifyGreaterThan(testCase, max(abs(ref.pred(21, :) - ra.pred(21, :))), ...
                  max(abs(rb.pred(21, :) - ra.pred(21, :))));
end

function testProtectedContractIsNeverFlagged(testCase)
[KS, x] = fixture(struct('spike', [20 3 0.003]));
protected = false(1, KS.nu); protected(3) = true;
r = innovationsLOO(KS, x, struct('protected', protected));
verifyFalse(testCase, any(r.flag(:, 3)));
end

function testMaxFlagFractionCapsTheDay(testCase)
[KS, x] = fixture(struct('spike', [20 1 0.003 ; 20 2 -0.003 ; 20 3 0.004]));
r = innovationsLOO(KS, x, struct('maxFlagFraction', 0.25));  % at most 2 of 8
verifyEqual(testCase, r.nFlag(20), 2);
verifyTrue(testCase, r.capped(20));
end

function testWrongParameterCount(testCase)
[KS, x] = fixture(struct());
verifyError(testCase, @() innovationsLOO(KS, x(1:end-1)), 'kalmanrates:ModelAlignment');
end

%% Fixture

function [KS, x] = fixture(opts)
% Level x_t (random walk) and nu contracts z_ti = g_i(x_t) + noise, with
% g_i(x) = (exp(x*T_i/365) - 1)/tc_i close to x (IBOR pricing in
% taylorInstrPC). opts.spike rows [day contract size] add to one quote,
% opts.jump [day size] moves the level.
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
tc = (T - 1)/365;  % EfInt(T) = (T-1)/365 with Ef = 1
g = @(lv) (exp(lv*(T-1)/365) - 1) ./ tc;

[KS.A, KS.B, KS.D, KS.G, KS.Iz, KS.z, KS.atAll, KS.oTAll, KS.oIndAll, KS.tcAll, KS.datesStepRelAll] = ...
  deal(cell(nt,1));
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
