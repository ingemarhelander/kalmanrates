function comparison = kalmanOptimizeSweep(parameterSpec, valuesOrMode)
% Runs kalmanOptimizeMain over one or more configuration parameters.
% Single parameter: kalmanOptimizeSweep('noisy_size', [1 3 10 100])
% Multiple parameters: kalmanOptimizeSweep(struct('noisy_size', [1 3], ...
%     'mu', [1e-3 1e-2]), 'cross')

isLegacyCall = ischar(parameterSpec) || ...
    (isstring(parameterSpec) && isscalar(parameterSpec));
if isLegacyCall
    parameterNames = {char(parameterSpec)};
    if ~isvarname(parameterNames{1})
        error('kalmanrates:InvalidSweepParameter', ...
            'parameterName must be a config key, such as ''noisy_size''.');
    end
    if nargin < 2
        error('kalmanrates:InvalidSweepValues', ...
            'Specify the values to use for "%s".', parameterNames{1});
    end
    rawValues = {valuesOrMode};
    mode = 'single';
elseif isstruct(parameterSpec) && isscalar(parameterSpec) && ...
        ~isempty(fieldnames(parameterSpec))
    parameterNames = fieldnames(parameterSpec)';
    rawValues = cellfun(@(name) parameterSpec.(name), parameterNames, ...
        'UniformOutput', false);
    mode = 'separate';
    if nargin >= 2
        mode = normalize_mode(valuesOrMode);
    end
else
    error('kalmanrates:InvalidSweepParameters', ...
        ['Specify one config key and its values, or a scalar struct mapping ' ...
         'config keys to value vectors.']);
end

scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
run(fullfile(projectRoot, 'setupPaths.m'));
config = readConfig(fullfile(scriptDir, 'kalmanOptimizeMain.config'));
parameterValues = cell(size(parameterNames));
defaults = cell(size(parameterNames));
for p = 1:numel(parameterNames)
    name = parameterNames{p};
    if ~isfield(config, name)
        error('kalmanrates:UnknownSweepParameter', ...
            'The config has no parameter named "%s".', name);
    end
    defaults{p} = config.(name);
    parameterValues{p} = normalize_values(rawValues{p}, defaults{p}, name);
end

runValues = build_run_values(parameterValues, defaults, mode);
numRuns = size(runValues, 1);
numParameters = numel(parameterNames);
objectives = nan(numRuns, 1);
gradientNorms = nan(numRuns, 1);
outOfSampleLogLik = nan(numRuns, 1);
outOfSamplePerDay = nan(numRuns, 1);
inSampleLogLik = nan(numRuns, 1);
histories = cell(numRuns, 1);
labels = cell(numRuns, 1);
mainScript = fullfile(scriptDir, 'kalmanOptimizeMain.m');

fprintf('Running %d optimization(s) in %s mode.\n', numRuns, mode);
for i = 1:numRuns
    labels{i} = run_label(parameterNames, runValues(i, :));
    fprintf('\n===== %s (%d of %d) =====\n', labels{i}, i, numRuns);
    override = cell2struct(runValues(i, :)', parameterNames', 1);
    saved = run_optimizer(mainScript, override);

    history = saved.history;
    histories{i} = struct('runs', history.runs, ...
        'bestRunIndex', history.bestRunIndex);
    bestRun = history.runs(history.bestRunIndex);
    objectives(i) = bestRun.fvals(end);
    gradientNorms(i) = bestRun.normOfgrads(end);
    if isfield(history, 'outOfSample')
        outOfSampleLogLik(i) = history.outOfSample.logLik;
        outOfSamplePerDay(i) = history.outOfSample.perDay;
        inSampleLogLik(i) = history.outOfSample.optimizeLogLik;
    end
end

metricNames = {'Objective', 'GradientNorm', 'OutOfSampleLogLik', 'InSampleLogLik'};
if isLegacyCall
    comparison = table(values_column(runValues(:, 1)), objectives, gradientNorms, ...
        outOfSampleLogLik, inSampleLogLik, ...
        'VariableNames', [{'ParameterValue'}, metricNames]);
else
    parameterColumns = cell(1, numParameters);
    for p = 1:numParameters
        parameterColumns{p} = values_column(runValues(:, p));
    end
    comparison = table(parameterColumns{:}, objectives, gradientNorms, ...
        outOfSampleLogLik, inSampleLogLik, ...
        'VariableNames', [parameterNames, metricNames]);
end

comparisonFile = unique_comparison_file(dataPath(), parameterNames, mode);
save(comparisonFile, 'comparison', 'parameterNames', 'parameterValues', ...
    'runValues', 'mode', 'histories', 'labels');
if isLegacyCall
    parameterName = parameterNames{1};
    values = valuesOrMode;
    save(comparisonFile, 'parameterName', 'values', '-append');
end
fprintf('\nComparison results:\n');
disp(comparison);
fprintf('Saved comparison to %s\n', comparisonFile);

plotData = struct('comparison', comparison, 'histories', {histories}, ...
    'labels', {labels}, 'outOfSamplePerDay', outOfSamplePerDay);
run_sweep_plots(fullfile(scriptDir, 'plots.m'), plotData);
end

function mode = normalize_mode(mode)
if isstring(mode) && isscalar(mode)
    mode = char(mode);
end
if ~ischar(mode) || ~isrow(mode)
    error('kalmanrates:InvalidSweepMode', ...
        'Mode must be "separate" or "cross".');
end
mode = lower(strtrim(mode));
if strcmp(mode, 'cartesian')
    mode = 'cross';
end
if ~any(strcmp(mode, {'separate', 'cross'}))
    error('kalmanrates:InvalidSweepMode', ...
        'Mode must be "separate" or "cross".');
end
end

function runValues = build_run_values(parameterValues, defaults, mode)
numParameters = numel(parameterValues);
if strcmp(mode, 'single')
    runValues = parameterValues{1};
elseif strcmp(mode, 'separate')
    numRuns = sum(cellfun(@numel, parameterValues));
    runValues = repmat(defaults, numRuns, 1);
    row = 0;
    for p = 1:numParameters
        for valueIndex = 1:numel(parameterValues{p})
            row = row + 1;
            runValues{row, p} = parameterValues{p}{valueIndex};
        end
    end
else
    indexValues = cellfun(@(values) 1:numel(values), parameterValues, ...
        'UniformOutput', false);
    indexGrids = cell(size(indexValues));
    [indexGrids{:}] = ndgrid(indexValues{:});
    numRuns = numel(indexGrids{1});
    runValues = cell(numRuns, numParameters);
    for p = 1:numParameters
        runValues(:, p) = parameterValues{p}(indexGrids{p}(:));
    end
end
end

function saved = run_optimizer(mainScript, override)
kalmanOptimizeMainConfigOverride = override; %#ok<NASGU>
run(mainScript);
saved = load(dataPath('kalmanOptimizeResult.mat'), 'history');
end

function run_sweep_plots(plotsScript, plotData)
kalmanOptimizeSweepPlotData = plotData; %#ok<NASGU>
run(plotsScript);
end

function sweepValues = normalize_values(values, configuredValue, parameterName)
if isnumeric(values) || islogical(values)
    if isempty(values) || ~isvector(values) || ~isreal(values) || ...
            any(~isfinite(double(values(:))))
        error('kalmanrates:InvalidSweepValues', ...
            'Numeric sweep values must be a nonempty finite vector.');
    end
    sweepValues = num2cell(values(:));
elseif isstring(values)
    if isempty(values) || ~isvector(values) || any(ismissing(values))
        error('kalmanrates:InvalidSweepValues', ...
            'Text sweep values must be a nonempty vector without missing values.');
    end
    sweepValues = cellstr(values(:));
elseif ischar(values) && isrow(values)
    sweepValues = {values};
elseif iscell(values) && ~isempty(values) && isvector(values)
    sweepValues = values(:);
else
    error('kalmanrates:InvalidSweepValues', ...
        'values must be a numeric, logical, or text vector.');
end

for i = 1:numel(sweepValues)
    value = sweepValues{i};
    if isstring(value) && isscalar(value)
        if ismissing(value)
            error('kalmanrates:InvalidSweepValue', ...
                'Sweep values for "%s" must not be missing.', parameterName);
        end
        value = char(value);
        sweepValues{i} = value;
    end
    if isnumeric(configuredValue)
        valid = isnumeric(value) && isscalar(value) && isreal(value) && isfinite(value);
    elseif islogical(configuredValue)
        valid = islogical(value) && isscalar(value);
    else
        valid = ischar(value) && isrow(value) && ~isempty(value);
    end
    if ~valid
        error('kalmanrates:InvalidSweepValue', ...
            'Sweep values for "%s" must match its configured value type.', parameterName);
    end
end
end

function column = values_column(sweepValues)
if all(cellfun(@(value) isnumeric(value) || islogical(value), sweepValues))
    column = cell2mat(sweepValues);
else
    column = sweepValues;
end
end

function label = run_label(parameterNames, runValues)
parts = cell(size(parameterNames));
for p = 1:numel(parameterNames)
    parts{p} = sprintf('%s = %s', parameterNames{p}, value_text(runValues{p}));
end
label = strjoin(parts, ', ');
end

function text = value_text(value)
if ischar(value)
    text = value;
else
    text = mat2str(value);
end
end

function file = unique_comparison_file(dataFolder, parameterNames, mode)
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
parameterLabel = strjoin(parameterNames, '_');
file = fullfile(dataFolder, sprintf( ...
    'kalmanOptimizeSweep_%s_%s_%s.mat', parameterLabel, mode, stamp));
suffix = 1;
while exist(file, 'file')
    file = fullfile(dataFolder, sprintf( ...
        'kalmanOptimizeSweep_%s_%s_%s_%d.mat', ...
        parameterLabel, mode, stamp, suffix));
    suffix = suffix + 1;
end
end
