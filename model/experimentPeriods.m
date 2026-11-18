function periods = experimentPeriods(config)
% experimentPeriods - In-sample, optimization and out-of-sample periods of an experiment.
%
%   periods = experimentPeriods(config)
%
%   config  struct from readConfig (kalmanOptimizeMain.config). Dates are
%           written yyyy-mm-dd and both end points are included. Keys:
%             insample_start, insample_end  days the model is built on
%                         (PCA and start parameters). Give both or neither.
%             optimize_start, optimize_end  days the parameters are
%                         optimized on. Must lie within in-sample; default
%                         whole in-sample.
%             out_of_sample_start, out_of_sample_end  out-of-sample period.
%                         out_of_sample_start must be after insample_end.
%                         Without out_of_sample_start there is no
%                         out-of-sample evaluation; without
%                         out_of_sample_end the last day of the data is used.
%             start_point init (the model's own start parameters) or saved
%                         (kalmanXoptAll.mat). Default init when in-sample
%                         is set, otherwise saved.
%
%   periods struct with the fields
%             hasInsample   true when insample_start/insample_end are set
%             insampleStart, insampleEnd, optimizeStart, optimizeEnd,
%             outOfSampleStart, outOfSampleEnd   datenums, [] when not given
%             buildEnd      last day the model is built for:
%                           outOfSampleEnd, or Inf (the last day of the data)
%             startPoint    'init' or 'saved'
%             modelFile     kalmanModel_<insample_start>_<insample_end>_<out_of_sample_end or end>.mat,
%                           or '' (the old kalmanDebug.mat) without in-sample
%
%   A broken rule gives the error kalmanrates:InvalidPeriods. Without
%   in-sample the old saved model is used and, when there is an
%   out-of-sample period, the warning kalmanrates:NotOutOfSample is given.

oldKeys = intersect(fieldnames(config), {'train_start', 'train_end', 'test_start', 'test_end'});
if ~isempty(oldKeys)
  fail(['%s: the keys are now called optimize_start/optimize_end and ' ...
        'out_of_sample_start/out_of_sample_end.'], strjoin(oldKeys, ', '));
end

periods = struct();
periods.insampleStart    = readDate(config, 'insample_start');
periods.insampleEnd      = readDate(config, 'insample_end');
periods.optimizeStart    = readDate(config, 'optimize_start');
periods.optimizeEnd      = readDate(config, 'optimize_end');
periods.outOfSampleStart = readDate(config, 'out_of_sample_start');
periods.outOfSampleEnd   = readDate(config, 'out_of_sample_end');

if isempty(periods.insampleStart) ~= isempty(periods.insampleEnd)
  fail('insample_start and insample_end must be given together.');
end
periods.hasInsample = ~isempty(periods.insampleStart);

checkOrder(periods.insampleStart, periods.insampleEnd, 'insample');
if ~isempty(periods.outOfSampleEnd) && isempty(periods.outOfSampleStart)
  fail('out_of_sample_end is given without out_of_sample_start.');
end
checkOrder(periods.outOfSampleStart, periods.outOfSampleEnd, 'out_of_sample');

if periods.hasInsample
  if isempty(periods.optimizeStart)
    periods.optimizeStart = periods.insampleStart;
  end
  if isempty(periods.optimizeEnd)
    periods.optimizeEnd = periods.insampleEnd;
  end
  if periods.optimizeStart < periods.insampleStart || periods.optimizeEnd > periods.insampleEnd
    fail('The optimization period %s to %s must lie within in-sample %s to %s.', ...
         ymd(periods.optimizeStart), ymd(periods.optimizeEnd), ...
         ymd(periods.insampleStart), ymd(periods.insampleEnd));
  end
  if ~isempty(periods.outOfSampleStart) && periods.outOfSampleStart <= periods.insampleEnd
    fail('out_of_sample_start %s must be after insample_end %s.', ...
         ymd(periods.outOfSampleStart), ymd(periods.insampleEnd));
  end
end
checkOrder(periods.optimizeStart, periods.optimizeEnd, 'optimize');

if isempty(periods.outOfSampleEnd)
  periods.buildEnd = Inf;
  endName = 'end';
else
  periods.buildEnd = periods.outOfSampleEnd;
  endName = ymd(periods.outOfSampleEnd);
end

if periods.hasInsample
  periods.modelFile = sprintf('kalmanModel_%s_%s_%s.mat', ymd(periods.insampleStart), ...
                              ymd(periods.insampleEnd), endName);
  defaultStart = 'init';
else
  periods.modelFile = '';
  defaultStart = 'saved';
  if ~isempty(periods.outOfSampleStart)
    warning('kalmanrates:NotOutOfSample', ...
            ['insample_start/insample_end are not set, so the old model (kalmanDebug.mat) ' ...
             'is used. Its PCA and start parameters may be built on the out-of-sample ' ...
             'days, so the evaluation is not guaranteed to be out of sample.']);
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
