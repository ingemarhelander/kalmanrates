function tests = testAdCounts
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
testMexDir = fullfile(fileparts(mfilename('fullpath')), 'mex');
addpath(testMexDir);
assert(exist('testAdCountsMex', 'file') == 3, ...
    'kalmanrates:MissingAdTestMex', ...
    'Build the test helper by running setupPaths and compileMex first.');
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testCountsAndDerivatives(testCase)
verifyEqual(testCase, testAdCountsMex(), 6);
end
