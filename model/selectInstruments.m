function [keep, labels] = selectInstruments(ricAll, nu, tenors)
% selectInstruments - Which of the model's instruments (contracts) to observe.
%
%   [keep, labels] = selectInstruments(ricAll, nu, tenors)
%
%   ricAll  RIC per instrument in model order (KI.ricAll). Only the first nu
%           are used; older saved models list instruments beyond 10Y that
%           were never part of the model.
%   nu      number of instruments in the model (K.nu).
%   tenors  'all', or a list (char with commas, string or cellstr) of
%           short names such as EONIA, 1W, 2W, 1M, 6M, 1Y, 18M, 10Y, or
%           full RICs such as EUREON1M=. Case does not matter.
%
%   keep    logical 1 x nu, true for the selected instruments.
%   labels  short name per instrument (1 x nu cellstr).

if numel(ricAll) < nu
  error('kalmanrates:ModelAlignment', ...
        'The model has %d instruments but only %d RICs are listed.', nu, numel(ricAll));
end
ricAll = cellstr(ricAll(1:nu));
labels = cellfun(@tenorLabel, ricAll, 'UniformOutput', false);

if nargin < 3 || isempty(tenors)
  tenors = 'all';
end
if ischar(tenors) || isstring(tenors)
  tenors = strsplit(char(tenors), ',');
end
tenors = upper(strtrim(cellstr(tenors)));
tenors = tenors(~cellfun(@isempty, tenors));

if isequal(tenors, {'ALL'})
  keep = true(1, nu);
  return;
end

keep = false(1, nu);
for i = 1:numel(tenors)
  match = strcmpi(labels, tenorLabel(tenors{i})) | strcmpi(ricAll, tenors{i});
  if ~any(match)
    error('kalmanrates:UnknownTenor', 'Unknown tenor "%s". Available: %s', ...
          tenors{i}, strjoin(labels, ', '));
  end
  keep = keep | match(:)';
end
end

function label = tenorLabel(name)
% EUREON3M= -> 3M, EUREONSW= -> 1W, EONIA= -> EONIA.
label = upper(strtrim(char(name)));
label = regexprep(label, '=$', '');
label = regexprep(label, '^EUREON', '');
if strcmp(label, 'SW')
  label = '1W';
end
end
