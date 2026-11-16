function [f, g] = kalmanObj(x, K, ~, ~)
% kalmanObj - Attrapp som ersätter model/kalmanObj.m i tester och demos.
%
% K är ett funktionshandtag [f, g] = K(x), så att line_search_main kan köras
% på en påhittad målfunktion utan Kalmanfilter, data eller mex.
% Läggs bara på sökvägen av de tester och demos som använder den.

    [f, g] = K(x);
end
