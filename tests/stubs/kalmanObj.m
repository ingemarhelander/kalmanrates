function [f, g] = kalmanObj(x, K, ~, ~)
% kalmanObj - TESTATTRAPP som ersätter model/kalmanObj.m i tester.
%
% K är ett funktionshandtag [f, g] = K(x), så att line_search_main kan köras
% på en påhittad målfunktion utan Kalmanfilter, data eller mex.
% Läggs bara på sökvägen i tester som behöver den (PathFixture).

    [f, g] = K(x);
end
