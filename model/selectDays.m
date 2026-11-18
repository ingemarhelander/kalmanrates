function indices = selectDays(times, startDate, endDate)
% selectDays - Indices of the observation days from startDate to endDate.
%
%   indices = selectDays(times, startDate, endDate)
%
%   times      datenum per observation (the time of day is ignored).
%   startDate  first day included (datenum), [] for the first day.
%   endDate    last day included (datenum), [] for the last day.
%
%   indices is a row vector of consecutive indices.

if isempty(startDate), startDate = -Inf; end
if isempty(endDate), endDate = Inf; end
days = floor(times(:));
indices = find(days >= startDate & days <= endDate)';
if isempty(indices)
  error('kalmanrates:NoObservations', 'No observation days between %s and %s.', ...
        showDate(startDate), showDate(endDate));
end
end

function s = showDate(d)
if isfinite(d)
  s = datestr(d, 'yyyy-mm-dd');
else
  s = 'the edge of the data';
end
end
