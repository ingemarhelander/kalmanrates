function indices = selectDays(times, opts)
% selectDays - Index of the observation days that match a date selection.
%
%   indices = selectDays(times, opts)
%
%   times  datenum per observation (the time of day is ignored).
%   opts   struct with the optional fields
%            startDate  first day included (datenum or 'yyyy-mm-dd'), default first day
%            endDate    last day included (datenum or 'yyyy-mm-dd'), default last day
%            weekdays   'all' or a list of days, e.g. 'wed' or 'mon,thu'
%                       (char, string or cellstr), default 'all'
%            dayStep    keep every dayStep:th day that remains after the
%                       date and weekday filters (1 = all), default 1
%
%   indices is a row vector in increasing time order.

if nargin < 2
  opts = struct();
end
days = floor(times(:));
if isempty(days)
  error('kalmanrates:NoObservations', 'No observations are available.');
end

startGiven = getOpt(opts, 'startDate', []);
endGiven = getOpt(opts, 'endDate', []);
startDate = parseDate(startGiven, days(1));
endDate = parseDate(endGiven, days(end));
if ~isempty(startGiven) && ~isempty(endGiven) && endDate < startDate
  error('kalmanrates:InvalidDateRange', 'endDate %s is before startDate %s.', ...
        datestr(endDate, 'yyyy-mm-dd'), datestr(startDate, 'yyyy-mm-dd'));
end
keep = days >= startDate & days <= endDate;

wanted = parseWeekdays(getOpt(opts, 'weekdays', 'all'));
if ~isempty(wanted)
  keep = keep & ismember(weekday(days), wanted);
end

indices = find(keep)';
dayStep = getOpt(opts, 'dayStep', 1);
if ~(isnumeric(dayStep) && isscalar(dayStep) && isfinite(dayStep) && ...
     dayStep >= 1 && dayStep == fix(dayStep))
  error('kalmanrates:InvalidDayStep', 'dayStep must be a positive integer.');
end
indices = indices(1:dayStep:end);

if isempty(indices)
  error('kalmanrates:NoObservations', ...
        'No observation days between %s and %s match the selection.', ...
        datestr(startDate, 'yyyy-mm-dd'), datestr(endDate, 'yyyy-mm-dd'));
end
end

function value = getOpt(opts, name, default)
value = default;
if isfield(opts, name) && ~isempty(opts.(name))
  value = opts.(name);
end
end

function d = parseDate(value, default)
if isempty(value)
  d = default;
elseif isnumeric(value) && isscalar(value)
  d = floor(value);
elseif ischar(value) || (isstring(value) && isscalar(value))
  try
    d = datenum(char(value), 'yyyy-mm-dd');
  catch
    error('kalmanrates:InvalidDate', 'Dates must be written yyyy-mm-dd, got "%s".', char(value));
  end
else
  error('kalmanrates:InvalidDate', 'Dates must be a datenum or a yyyy-mm-dd string.');
end
end

function wanted = parseWeekdays(value)
% Returns MATLAB weekday numbers (Sun = 1), or [] for all days.
names = {'sun', 'mon', 'tue', 'wed', 'thu', 'fri', 'sat'};
if ischar(value) || isstring(value)
  value = strsplit(char(value), ',');
end
value = lower(strtrim(cellstr(value)));
if isequal(value, {'all'})
  wanted = [];
  return;
end
wanted = zeros(1, numel(value));
for i = 1:numel(value)
  k = find(strncmp(names, value{i}, 3) & strlength(value{i}) >= 3, 1);
  if isempty(k)
    error('kalmanrates:InvalidWeekday', ...
          'Unknown weekday "%s" (use all or mon, tue, wed, thu, fri).', value{i});
  end
  wanted(i) = k;
end
end
