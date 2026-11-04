% Build the likelihood and its AD regression test from any working directory.
% Requires a configured MATLAB C++ compiler and the supervisor's libraries.
mexDir = fileparts(mfilename('fullpath'));
testMexDir = fullfile(mexDir, '..', 'tests', 'mex');

% Same library locations as before, now resolved relative to this script.
% Adjust progDir if the libraries are installed elsewhere on this computer.
progDir = fullfile(mexDir, '..', '..', '..', '..', 'prog');
boostDir = fullfile(progDir, 'diverse', 'boost_1_67_0');
mexOptions = {['-I' mexDir], ['-I' boostDir], ...
    ['-I' fullfile(progDir, 'XiMathBoost', 'include')], ...
    ['-I' fullfile(progDir, 'XiFinBoost', 'include')], ...
    ['-I' fullfile(progDir, 'XiFinPortfolio', 'include')], ...
    ['-L' fullfile(boostDir, 'stage', 'lib')], ...
    '-lboost_date_time-vc141-mt-x64-1_67'};

% Unload old binaries before replacing them, particularly on Windows.
clear mexKalmanLogLikelihood testAdCountsMex
mex(mexOptions{:}, '-outdir', mexDir, ...
    fullfile(mexDir, 'mexKalmanLogLikelihood.cpp'));
mex(mexOptions{:}, '-outdir', testMexDir, ...
    fullfile(testMexDir, 'testAdCountsMex.cpp'));
rehash
fprintf('Built the likelihood in %s and the AD test in %s.\n', mexDir, testMexDir);
which mexKalmanLogLikelihood -all

% Historical compiler examples below are not executed.
return

% mex  -I../../../../prog/boost_1_57_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/boost_1_57_0/stage/lib mexPortfolio.cpp
%mex  -I../../../../prog/boost_1_57_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/boost_1_57_0/bin.v2/libs/date_time/build/msvc-10.0/release/link-static/threading-multi -lboost_date_time-vc100-mt-1_57 mexPortfolio.cpp
% mex  -I../../../../prog/diverse/boost_1_63_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinPortfolio/include mexPortfolio.cpp
% mex  -I../../../../prog/diverse/boost_1_49_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_49_0/bin.v2/libs/date_time/build/msvc-10.0/release/link-static/threading-multi -lboost_date_time-vc100-mt-1_57 mexPortfolio.cpp

% Windows:
% mex  -g -I../../../../prog/diverse/boost_1_67_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_67_0/stage/lib -lboost_date_time-vc141-mt-x64-1_67 mexPortfolio.cpp
mex  -I../../../../prog/diverse/boost_1_67_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_67_0/stage/lib -lboost_date_time-vc141-mt-x64-1_67 mexPortfolio.cpp

% Linux:
% mex  -g -I../../../../prog/diverse/boost_1_67_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_67_0/stage/lib mexPortfolio.cpp
% mex  -I../../../../prog/diverse/boost_1_67_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_67_0/stage/lib mexPortfolio.cpp


