% Run from the mex/ folder (cd mex) - include/library paths are relative to it.
mex  -I../../../../prog/diverse/boost_1_67_0 -I../../../../prog/XiMathBoost/include -I../../../../prog/XiFinBoost/include -I../../../../prog/XiFinPortfolio/include -L../../../../prog/diverse/boost_1_67_0/stage/lib -lboost_date_time-vc141-mt-x64-1_67 mexKalmanLogLikelihood.cpp
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


