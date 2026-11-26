function name = useDataSet(target)
% useDataSet - Choose which market data is in data/EUR.mat and data/EUR100.mat.
%
%   useDataSet('raw')     the original data (from createMatFromExcel)
%   useDataSet('clean')   the data where cleanMarketData removed noisy quotes
%   useDataSet('backup')  save the active (raw) data as EUR_raw.mat and
%                         EUR100_raw.mat if they do not exist yet
%   name = useDataSet()   print and return the active data set
%
%   The existing scripts always read EUR.mat and EUR100.mat, so the data sets
%   are kept as copies next to them and copied in when chosen:
%
%     data/EUR_raw.mat     data/EUR100_raw.mat     original
%     data/EUR_clean.mat   data/EUR100_clean.mat   cleaned (EUR100_clean.mat
%                                                  exists once
%                                                  runTenorMultiIrCurves has
%                                                  run on the cleaned data)
%
%   Without EUR100_clean.mat, 'clean' uses the raw curves (EUR100_raw.mat):
%   only the quotes of the Kalman model are cleaned, the PCA is unchanged.
%
%   When runTenorMultiIrCurves has written a new EUR100.mat while 'clean' is
%   active, it is saved as EUR100_clean.mat before switching. A file that
%   has changed in any other way is never overwritten; the function stops
%   and says which file.
%
%   The active set is recorded in data/dataSet.mat.

  if nargin < 1
    st = readState();
    fprintf('Active data set: %s\n', st.name);
    name = st.name;
    return;
  end
  target = lower(char(target));

  switch target
    case 'backup'
      backup();
    case {'raw', 'clean'}
      switchTo(target);
    otherwise
      error('kalmanrates:InvalidDataSet', 'Use useDataSet(''raw''), (''clean'') or (''backup'').');
  end
  if nargout > 0
    name = readState().name;
  end
end

function backup()
  st = readState();
  if ~strcmp(st.name, 'raw')
    error('kalmanrates:InvalidDataSet', 'The active data set is %s, run useDataSet(''raw'') first.', st.name);
  end
  for base = {'EUR', 'EUR100'}
    store = dataPath([base{1} '_raw.mat']);
    if isfile(store)
      if hasChanged(st, base{1})
        error('kalmanrates:DataSetChanged', ...
              ['%s.mat differs from %s, which already exists. Remove or rename one of ' ...
               'them by hand so that the raw data set is unambiguous.'], base{1}, store);
      end
    else
      fprintf('Saving %s.mat as %s\n', base{1}, store);
      copyfile(dataPath([base{1} '.mat']), store);
    end
  end
  writeState('raw');
end

function switchTo(target)
  st = readState();

  % Keep what has been produced in the active set
  if hasChanged(st, 'EUR')
    if strcmp(st.name, 'raw') && ~isfile(dataPath('EUR_raw.mat'))
      fprintf('Saving EUR.mat as %s\n', dataPath('EUR_raw.mat'));
      copyfile(dataPath('EUR.mat'), dataPath('EUR_raw.mat'));
    else
      error('kalmanrates:DataSetChanged', ...
            ['EUR.mat has changed since the %s data set was chosen. If it was built again ' ...
             'with createMatFromExcel, remove EUR_raw.mat and EUR_clean.mat by hand and ' ...
             'run cleanMarketData again.'], st.name);
    end
  end
  if hasChanged(st, 'EUR100')
    switch st.name
      case 'clean'
        fprintf('Saving EUR100.mat (curves of the cleaned data) as %s\n', dataPath('EUR100_clean.mat'));
        copyfile(dataPath('EUR100.mat'), dataPath('EUR100_clean.mat'));
      case 'raw'
        if isfile(dataPath('EUR100_raw.mat'))
          error('kalmanrates:DataSetChanged', ...
                ['EUR100.mat has changed since the raw data set was chosen and EUR100_raw.mat ' ...
                 'already exists. Rename one of them by hand.']);
        end
        fprintf('Saving EUR100.mat as %s\n', dataPath('EUR100_raw.mat'));
        copyfile(dataPath('EUR100.mat'), dataPath('EUR100_raw.mat'));
    end
  end

  % Copy in the chosen set
  eurSource = dataPath(['EUR_' target '.mat']);
  if ~isfile(eurSource)
    error('kalmanrates:MissingDataSet', '%s does not exist. Run cleanMarketData first.', eurSource);
  end
  curveSource = dataPath(['EUR100_' target '.mat']);
  if ~isfile(curveSource)
    if strcmp(target, 'clean')
      curveSource = dataPath('EUR100_raw.mat');
      fprintf(['No EUR100_clean.mat: the raw curves are used (PCA unchanged). Run ' ...
               'runTenorMultiIrCurves to build curves of the cleaned data.\n']);
    end
    if ~isfile(curveSource)
      error('kalmanrates:MissingDataSet', '%s does not exist.', curveSource);
    end
  end
  copyIfDifferent(eurSource, dataPath('EUR.mat'));
  copyIfDifferent(curveSource, dataPath('EUR100.mat'));
  writeState(target);
  [~, e] = fileparts(eurSource); [~, c] = fileparts(curveSource);
  fprintf('Active data set: %s (EUR.mat = %s.mat, EUR100.mat = %s.mat)\n', target, e, c);
end

function copyIfDifferent(source, dest)
  if isfile(dest) && sameFile(dir(source), dir(dest))
    return;
  end
  copyfile(source, dest);
end

function same = sameFile(a, b)
  same = a.bytes == b.bytes && abs(a.datenum - b.datenum) < 1/86400;
end

function changed = hasChanged(st, base)
% True when data/<base>.mat is not the file recorded when the set was chosen.
  file = dataPath([base '.mat']);
  if ~isfile(file)
    changed = false;
    return;
  end
  if ~isfield(st.files, base)
    changed = true;  % never recorded
    return;
  end
  d = dir(file);
  changed = ~(d.bytes == st.files.(base).bytes && abs(d.datenum - st.files.(base).datenum) < 1/86400);
end

function st = readState()
  file = dataPath('dataSet.mat');
  if isfile(file)
    st = load(file);
  else
    st = struct('name', 'raw', 'files', struct());
  end
end

function writeState(name)
  files = struct();
  for base = {'EUR', 'EUR100'}
    file = dataPath([base{1} '.mat']);
    if isfile(file)
      d = dir(file);
      files.(base{1}) = struct('bytes', d.bytes, 'datenum', d.datenum);
    end
  end
  save(dataPath('dataSet.mat'), 'name', 'files');
end
