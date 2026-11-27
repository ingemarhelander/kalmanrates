function tests = testKalmanOptimizeSweep
tests = functiontests(localfunctions);
end

function setupOnce(~)
run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
end

function testInvalidModeIsRejected(testCase)
parameters = struct('noisy_size', [1 3], 'mu', [0.01 0.1]);
verifyError(testCase, @() kalmanOptimizeSweep(parameters, 'diagonal'), ...
    'kalmanrates:InvalidSweepMode');
end

function testPlotCreatesComparisonFigures(testCase)
oldVisibility = get(groot, 'defaultFigureVisible');
figuresBefore = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() cleanupFigures(figuresBefore, oldVisibility)); %#ok<NASGU>
set(groot, 'defaultFigureVisible', 'off');

runData = struct('fvals', [1; 2; 3], 'normOfgrads', [4; 2; 1]);
plotData.histories = {struct('runs', runData, 'bestRunIndex', 1)};
plotData.labels = {'noisy_size = 1'};
plotData.comparison = table(3, 2, 'VariableNames', {'Objective', 'AverageDSize'});
plotData.outOfSamplePerDay = 0.2;
kalmanOptimizeSweepPlotData = plotData; %#ok<NASGU>
run(fullfile(fileparts(mfilename('fullpath')), '..', 'scripts', 'plots.m'));

newFigures = setdiff(findall(groot, 'Type', 'figure'), figuresBefore);
verifyEqual(testCase, numel(newFigures), 2);
end

function cleanupFigures(figuresBefore, oldVisibility)
newFigures = setdiff(findall(groot, 'Type', 'figure'), figuresBefore);
delete(newFigures);
set(groot, 'defaultFigureVisible', oldVisibility);
end
