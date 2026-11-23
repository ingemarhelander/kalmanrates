function config = readConfig(configPath)
% readConfig - Read key/value rows of a config file into a struct.
%
%   config = readConfig(configPath)
%
%   Each row is "key value". Blank lines and lines beginning with # or %
%   are ignored. true/false become logicals, numbers become doubles and
%   everything else (for example dates yyyy-mm-dd) is kept as text.

fid = fopen(configPath, 'r');
if fid == -1
  error('kalmanrates:MissingConfig', 'Could not open config file: %s', configPath);
end
cleanup = onCleanup(@() fclose(fid));
config = struct();

while true
  line = fgetl(fid);
  if ~ischar(line)
    break;
  end
  line = strtrim(line);
  if isempty(line) || startsWith(line, '#') || startsWith(line, '%')
    continue;
  end

  tokens = regexp(line, '^([A-Za-z]\w*)\s+(.+)$', 'tokens', 'once');
  if isempty(tokens) || ~isvarname(tokens{1})
    error('kalmanrates:InvalidConfig', 'Invalid config row: %s', line);
  end

  key = tokens{1};
  valueText = strtrim(tokens{2});
  if strcmpi(valueText, 'true')
    value = true;
  elseif strcmpi(valueText, 'false')
    value = false;
  else
    numericValue = str2double(valueText);
    if ~isnan(numericValue) || strcmpi(valueText, 'nan')
      value = numericValue;
    else
      value = valueText;
    end
  end
  config.(key) = value;
end
end
