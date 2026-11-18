function periods = experimentPeriods(config)
% experimentPeriods - In-sample, training and test periods of an experiment.
%
%   periods = experimentPeriods(config)
%
%   config  struct from readConfig (kalmanOptimizeMain.config). Dates are
%           written yyyy-mm-dd and both end points are included. Keys:
%             insample_start, insample_end  days the model is built on (PCA
%                         and start parameters). Give both or neither.
%             train_start, train_end  days the parameters are optimized on.
%                         Must lie within in-sample; default whole in-sample.
%             test_start, test_end  out-of-sample period. test_start must be
%                         after insample_end. Without test_start there is no
%                         test; without test_end the last day of the data.
%             start_point init (the model's own start parameters) or saved
%                         (kalmanXoptAll.mat). Default init when in-sample
%                         is set, otherwise saved.
%
%   periods struct with the fields
%             hasInsample   true when insample_start/insample_end are set
%             insampleStart, insampleEnd, trainStart, trainEnd,
%             testStart, testEnd   datenums, [] when not given
%             buildEnd      last day the model is built for: testEnd, or
%                           Inf (the last day of the data)
%             startPoint    'init' or 'saved'
%             modelFile     kalmanModel_<insample_start>_<insample_end>_<test_end or end>.mat,
%                           or '' (the old kalmanDebug.mat) without in-sample
%
%   A broken rule gives the error kalmanrates:InvalidPeriods. Without
%   in-sample the old saved model is used and, when there is a test period,
%   the warning kalmanrates:NotOutOfSample is given.

periods = struct();
periods.insampleStart = readDate(config, 'insample_start');
periods.insampleEnd   = readDate(config, 'insample_end');
periods.trainStart    = readDate(config, 'train_start');
periods.trainEnd      = readDate(config, 'train_end');
periods.testStart     = readDate(config, 'test_start');
periods.testEnd       = readDate(config, 'test_end');

if isempty(periods.insampleStart) ~= isempty(periods.insampleEnd)
  fail('insample_start and insample_end must be given together.');
end
periods.hasInsample = ~isempty(periods.insampleStart);

checkOrder(periods.insampleStart, periods.insampleEnd, 'insample');
if ~isempty(periods.testEnd) && isempty(periods.testStart)
  fail('test_end is given without test_start.');
end
checkOrder(periods.testStart, periods.testEnd, 'test');

if periods.hasInsample
  if isempty(periods.trainStart)
    periods.trainStart = periods.insampleStart;
  end
  if isempty(periods.trainEnd)
    periods.trainEnd = periods.insampleEnd;
  end
  if periods.trainStart < periods.insampleStart || periods.trainEnd > periods.insampleEnd
    fail('The training period %s to %s must lie within in-sample %s to %s.', ...
         ymd(periods.trainStart), ymd(periods.trainEnd), ...
         ymd(periods.insampleStart), ymd(periods.insampleEnd));
  end
  if ~isempty(periods.testStart) && periods.testStart <= periods.insampleEnd
    fail('test_start %s must be after insample_end %s.', ...
         ymd(periods.testStart), ymd(periods.insampleEnd));
  end
end
checkOrder(periods.trainStart, periods.trainEnd, 'train');

if isempty(periods.testEnd)
  periods.buildEnd = Inf;
  endName = 'end';
else
  periods.buildEnd = periods.testEnd;
  endName = ymd(periods.testEnd);
end

if periods.hasInsample
  periods.modelFile = sprintf('kalmanModel_%s_%s_%s.mat', ymd(periods.insampleStart), ...
                              ymd(periods.insampleEnd), endName);
  defaultStart = 'init';
else
  periods.modelFile = '';
  defaultStart = 'saved';
  if ~isempty(periods.testStart)
    warning('kalmanrates:NotOutOfSample', ...
            ['insample_start/insample_end are not set, so the old model (kalmanDebug.mat) ' ...
             'is used. Its PCA and start parameters may be built on the test days, so ' ...
             'the test is not guaranteed to be out of sample.']);
  end
end

periods.startPoint = defaultStart;
if isfield(config, 'start_point')
  periods.startPoint = lower(strtrim(char(string(config.start_point))));
end
if ~ismember(periods.startPoint, {'init', 'saved'})
  fail('start_point must be init or saved, got "%s".', periods.startPoint);
end
end

function d = readDate(config, key)
d = [];
if ~isfield(config, key) || isempty(config.(key))
  return;
end
text = config.(key);
if ~(ischar(text) || (isstring(text) && isscalar(text))) || ...
   isempty(regexp(char(text), '^\d{4}-\d{2}-\d{2}$', 'once'))
  fail('%s must be a date written yyyy-mm-dd.', key);
end
text = char(text);
d = datenum(text, 'yyyy-mm-dd');
if ~strcmp(ymd(d), text)
  fail('%s is not a valid date: %s.', key, text);
end
end

function checkOrder(first, last, name)
if ~isempty(first) && ~isempty(last) && last < first
  fail('%s_end %s is before %s_start %s.', name, ymd(last), name, ymd(first));
end
end

function s = ymd(d)
s = datestr(d, 'yyyy-mm-dd');
end

function fail(varargin)
error('kalmanrates:InvalidPeriods', varargin{:});
end
