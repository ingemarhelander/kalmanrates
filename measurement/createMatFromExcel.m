% function createMatFromExcel(fileName)
% fileName = 'CHF.xlsx';
fileName = 'EUR.xlsx';
% fileName = 'EUR_ICAP.xlsx';
% fileName = 'EUR_EONIA.xlsx';
% fileName = 'GBP.xlsx';
% fileName = 'JPY.xlsx';
% fileName = 'KRW.xlsx';
% fileName = 'SEK.xlsx';
% fileName = 'USD.xlsx';
% fileName = 'USDSEK.xlsx';
% fileName = 'USDNOK.xlsx';
% fileName = 'EURSEK.xlsx';
% fileName = 'EURUSD.xlsx';
% fileName = 'GBPUSD.xlsx';
% fileName = 'AUDUSD.xlsx';
% fileName = 'USDCAD.xlsx';
% fileName = 'USDCHF.xlsx';
% fileName = 'USDJPY.xlsx';
% fileName = 'USDZAR.xlsx';
% fileName = 'EUA.xlsx';

[type, sheetNames]=xlsfinfo(fileName);

i = 1; % Counter to deal with AFOSHEET
clear data;

for k = 1:length(sheetNames)
  if (strcmp(sheetNames(k), 'settings'))
    continue;
  end
%   [data{i},txt{i},raw{i}] = xlsread(fileName, char(sheetNames(k)));
  d = readtable(fileName, "Sheet", char(sheetNames(k)));

  sheet{i} = sheetNames(k);
  ric{i} = d.Properties.VariableDescriptions(2:end);

  dates{i} = datenum(table2array(d(:,1)));
  c = table2cell(d(:,2:end));
  data{i} = ones(size(c))*NaN;
  for j=1:size(c,2)
    ind = cellfun(@(x) ischar(x), c(:,j));
    c(ind,j) = cell(1,1);
    tmp = cell2mat(c(:,j)); % cell2mat removes empty columns, ensures that columns are kept
    if (~isempty(tmp))
      data{i}(:,j) = tmp;
    end
  end
%   data{i} = cell2mat(c);
  if (length(dates{i}) > 1)
    if (dates{i}(2) < dates{i}(1)) % Dates are sorted in decreasing order, change to increasing order
      n = size(data{i},1);
      indData = n:-1:1;
      data{i} = data{i}(indData,:);
      dates{i} = dates{i}(indData,1);
    end
  end
  i = i+1;
end

k = strfind(fileName, '.');
fileNameMatlab = fileName(1:k(end)-1);

save(fileNameMatlab, 'sheet', 'ric', 'dates', 'data', 'sheetNames');

