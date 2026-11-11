function compileMex(boostDir)
% Build the likelihood and AD test with local sparseMatrix.h and Boost headers.
%   compileMex('C:/path/to/boost_1_67_0')
% Without an argument, use BOOST_ROOT or ~/Documents/boost_1_67_0 on Windows.
mexDir = fileparts(mfilename('fullpath'));
testMexDir = fullfile(fileparts(mexDir), 'tests', 'mex');
if nargin < 1 || isempty(boostDir)
  boostDir = getenv('BOOST_ROOT');
  if isempty(boostDir)
    boostDir = fullfile(getenv('USERPROFILE'), 'Documents', 'boost_1_67_0');
  end
end
boostDir = char(boostDir);
assert(isfile(fullfile(boostDir, 'boost', 'numeric', 'ublas', 'matrix.hpp')) && ...
       isfile(fullfile(boostDir, 'boost', 'lexical_cast.hpp')), ...
    'kalmanrates:MissingBoost', ...
    'Pass the unpacked Boost directory to compileMex, or set BOOST_ROOT.');
assert(~isempty(mex.getCompilerConfigurations('C++', 'Selected')), ...
    'kalmanrates:MissingCppCompiler', 'Configure a compiler using mex -setup cpp.');

% uBLAS and lexical_cast are header-only. No XiMath/XiFin or date_time library
% is linked; the sparse matrix implementation comes from this repository.
mexOptions = {['-I' mexDir], ['-I' boostDir]};
clear mexKalmanLogLikelihood testAdCountsMex
% Build the small helper first, before replacing the production binary.
mex(mexOptions{:}, '-outdir', testMexDir, ...
    fullfile(testMexDir, 'testAdCountsMex.cpp'));
mex(mexOptions{:}, '-outdir', mexDir, ...
    fullfile(mexDir, 'mexKalmanLogLikelihood.cpp'));
rehash
fprintf('Built likelihood and AD test using Boost in %s.\n', boostDir);
which mexKalmanLogLikelihood -all
end
