% setupPaths - Lägger till projektets mappar på MATLAB-sökvägen.
% Anropas automatiskt i början av skripten i scripts/, men kan även köras
% manuellt en gång per MATLAB-session.

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);
addpath(fullfile(projectRoot, 'scripts'));
addpath(fullfile(projectRoot, 'model'));
addpath(fullfile(projectRoot, 'mex'));
addpath(fullfile(projectRoot, 'measurement'));
addpath(fullfile(projectRoot, 'data'));
addpath(fullfile(projectRoot, 'data', 'raw'));
clear projectRoot
