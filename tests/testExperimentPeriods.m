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
                          'test_start', '2019-01-01', 'test_end', '2019-03-31'));
verifyTrue(testCase, p.hasInsample);
verifyEqual(testCase, p.trainStart, datenum(2016,1,1));
verifyEqual(testCase, p.trainEnd, datenum(2018,12,31));
verifyEqual(testCase, p.testStart, datenum(2019,1,1));
verifyEqual(testCase, p.buildEnd, datenum(2019,3,31));
verifyEqual(testCase, p.startPoint, 'init');
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_2019-03-31.mat');
end

function testNoTestEndUsesEnd(testCase)
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'test_start', '2019-01-01'));
verifyEmpty(testCase, p.testEnd);
verifyEqual(testCase, p.buildEnd, Inf);
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_end.mat');
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31'));
verifyEmpty(testCase, p.testStart);
verifyEqual(testCase, p.modelFile, 'kalmanModel_2016-01-01_2018-12-31_end.mat');
end

function testTrainWithinInsample(testCase)
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'train_start', '2018-01-01', 'test_start', '2019-01-01', ...
                          'start_point', 'saved'));
verifyEqual(testCase, p.trainStart, datenum(2018,1,1));
verifyEqual(testCase, p.trainEnd, datenum(2018,12,31));
verifyEqual(testCase, p.startPoint, 'saved');
% The whole in-sample, end points included, is allowed.
p = experimentPeriods(cfg('insample_start', '2016-01-01', 'insample_end', '2018-12-31', ...
                          'train_start', '2016-01-01', 'train_end', '2018-12-31'));
verifyEqual(testCase, [p.trainStart p.trainEnd], datenum([2016 2018], [1 12], [1 31]));
end

function testWithoutInsampleWorksAsBefore(testCase)
p = verifyWarning(testCase, @() experimentPeriods(cfg('train_start', '2018-01-01', ...
        'train_end', '2018-12-31', 'test_start', '2019-01-01')), 'kalmanrates:NotOutOfSample');
verifyFalse(testCase, p.hasInsample);
verifyEqual(testCase, p.modelFile, '');
verifyEqual(testCase, p.startPoint, 'saved');
verifyEqual(testCase, p.trainStart, datenum(2018,1,1));
% Without a test there is nothing to warn about, and missing dates stay empty.
p = verifyWarningFree(testCase, @() experimentPeriods(struct()));
verifyEmpty(testCase, p.trainStart);
verifyEmpty(testCase, p.trainEnd);
verifyEmpty(testCase, p.testStart);
% Without in-sample the test may overlap the training (evaluateOutOfSample warns).
verifyWarning(testCase, @() experimentPeriods(cfg('train_start', '2018-01-01', ...
        'test_start', '2017-01-01')), 'kalmanrates:NotOutOfSample');
end

%% Broken rules

function testInvalidPeriods(testCase)
ins = {'insample_start', '2016-01-01', 'insample_end', '2018-12-31'};
bad = {
    {'insample_start', '2016-01-01'}                          % end missing
    {'insample_end', '2018-12-31'}                            % start missing
    {'insample_start', '2018-01-01', 'insample_end', '2017-12-31'}
    [ins, {'train_start', '2015-12-31'}]                      % before in-sample
    [ins, {'train_end', '2019-01-01'}]                        % after in-sample
    [ins, {'train_start', '2018-06-01', 'train_end', '2018-01-01'}]
    [ins, {'test_start', '2018-12-31'}]                       % test in in-sample
    [ins, {'test_start', '2018-06-01'}]
    [ins, {'test_start', '2019-02-01', 'test_end', '2019-01-31'}]
    [ins, {'test_end', '2019-03-31'}]                         % no test_start
    [ins, {'start_point', 'xyz'}]
    [ins, {'insample_start', '2016/01/01'}]                   % date format
    [ins, {'train_start', '2018-02-30'}]                      % no such day
    [ins, {'train_end', 20181231}]                            % number, not date
    {'train_start', '2018-06-01', 'train_end', '2018-01-01'}  % without in-sample
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
