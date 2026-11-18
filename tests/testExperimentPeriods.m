function tests = testExperimentPeriods
% Tests for experimentPeriods and readConfig.
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

%% Defaults and file name

function testInsampleDefaults(testCase)
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'out_of_sample_start', '2019-01-01', 'out_of_sample_end', '2019-03-31'));
verifyTrue(testCase, p.hasInsample);
verifyEqual(testCase, p.optimizeStart, datenum(2016,1,1));
verifyEqual(testCase, p.optimizeEnd, datenum(2018,12,31));
verifyEqual(testCase, p.outOfSampleStart, datenum(2019,1,1));
verifyEqual(testCase, p.buildEnd, datenum(2019,3,31));
verifyEqual(testCase, p.startPoint, 'init');
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_2019-03-31.mat');
end

function testNoOutOfSampleEndUsesEnd(testCase)
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'out_of_sample_start', '2019-01-01'));
verifyEmpty(testCase, p.outOfSampleEnd);
verifyEqual(testCase, p.buildEnd, Inf);
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_end.mat');
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31'));
verifyEmpty(testCase, p.outOfSampleStart);
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_end.mat');
end

function testOptimizeWithinInsample(testCase)
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'optimize_start', '2018-01-01', 'out_of_sample_start', '2019-01-01', ...
                          'start_point', 'saved'));
verifyEqual(testCase, p.optimizeStart, datenum(2018,1,1));
verifyEqual(testCase, p.optimizeEnd, datenum(2018,12,31));
verifyEqual(testCase, p.startPoint, 'saved');
% The whole in-sample, end points included, is allowed.
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'optimize_start', '2016-01-01', 'optimize_end', '2018-12-31'));
verifyEqual(testCase, [p.optimizeStart p.optimizeEnd], datenum([2016 2018], [1 12], [1 31]));
end

function testWithoutInsampleWorksAsBefore(testCase)
p = verifyWarning(testCase, @() experimentPeriods(cfg('optimize_start', '2018-01-01', ...
        'optimize_end', '2018-12-31', 'out_of_sample_start', '2019-01-01')), 'kalmanrates:NotOutOfSample');
verifyFalse(testCase, p.hasInsample);
verifyEqual(testCase, p.modelFile, '');
verifyEqual(testCase, p.startPoint, 'saved');
verifyEqual(testCase, p.optimizeStart, datenum(2018,1,1));
% Without an out-of-sample period there is nothing to warn about, and
% missing dates stay empty.
p = verifyWarningFree(testCase, @() experimentPeriods(struct()));
verifyEmpty(testCase, p.optimizeStart);
verifyEmpty(testCase, p.optimizeEnd);
verifyEmpty(testCase, p.outOfSampleStart);
% Without in-sample the out-of-sample period may overlap the optimization
% period (evaluateOutOfSample then stops with an error).
verifyWarning(testCase, @() experimentPeriods(cfg('optimize_start', '2018-01-01', ...
        'out_of_sample_start', '2017-01-01')), 'kalmanrates:NotOutOfSample');
end

%% Broken rules

function testInvalidPeriods(testCase)
ins = {'insample_start', '2016-01-01', 'insample_end', '2018-12-31'};
bad = {
    {'insample_start', '2016-01-01'}                          % end missing
    {'insample_end', '2018-12-31'}                            % start missing
    {'insample_start', '2018-01-01', 'insample_end', '2017-12-31'}
    [ins, {'optimize_start', '2015-12-31'}]                   % before in-sample
    [ins, {'optimize_end', '2019-01-01'}]                     % after in-sample
    [ins, {'optimize_start', '2018-06-01', 'optimize_end', '2018-01-01'}]
    [ins, {'out_of_sample_start', '2018-12-31'}]              % out-of-sample in in-sample
    [ins, {'out_of_sample_start', '2018-06-01'}]
    [ins, {'out_of_sample_start', '2019-02-01', 'out_of_sample_end', '2019-01-31'}]
    [ins, {'out_of_sample_end', '2019-03-31'}]                % no out_of_sample_start
    [ins, {'start_point', 'xyz'}]
    [ins, {'insample_start', '2016/01/01'}]                   % date format
    [ins, {'optimize_start', '2018-02-30'}]                   % no such day
    [ins, {'optimize_end', 20181231}]                         % number, not date
    {'optimize_start', '2018-06-01', 'optimize_end', '2018-01-01'}  % without in-sample
    [ins, {'train_start', '2018-01-01'}]                      % old key names
    [ins, {'test_start', '2019-01-01'}]
};
for i = 1:numel(bad)
  verifyError(testCase, @() experimentPeriods(cfg(bad{i}{:})), 'kalmanrates:InvalidPeriods', ...
              sprintf('Case %d should be rejected.', i));
end
end

%% readConfig

function testReadConfigKeepsDatesAsText(testCase)
file = [tempname '.config'];
cleanup = onCleanup(@() delete(file));
fid = fopen(file, 'w');
fprintf(fid, '# comment\n%% comment\n\ninsample_start 2016-01-01\nmaxIter 100\nmu 1E-2\nflag true\nlabel some text\n');
fclose(fid);
c = readConfig(file);
verifyEqual(testCase, c.insample_start, '2016-01-01');
verifyEqual(testCase, c.maxIter, 100);
verifyEqual(testCase, c.mu, 1e-2);
verifyEqual(testCase, c.flag, true);
verifyEqual(testCase, c.label, 'some text');
% insample_end is missing in this file.
verifyError(testCase, @() experimentPeriods(c), 'kalmanrates:InvalidPeriods');
end

function testShippedConfigIsValid(testCase)
file = fullfile(fileparts(mfilename('fullpath')), '..', 'scripts', 'kalmanOptimizeMain.config');
c = readConfig(file);
verifyWarningFree(testCase, @() experimentPeriods(c));
end

%% Helpers

function c = cfg(varargin)
% Config struct from key/value pairs; a later key overrides an earlier one.
c = struct();
for i = 1:2:numel(varargin)
  c.(varargin{i}) = varargin{i+1};
end
end
