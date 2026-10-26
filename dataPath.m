function p = dataPath(varargin)
% dataPath - Fullständig sökväg till en fil i projektets data-mapp.
%   dataPath('EUR100.mat')          -> <projekt>/data/EUR100.mat
%   dataPath('raw', 'EUR.xlsx')     -> <projekt>/data/raw/EUR.xlsx

p = fullfile(fileparts(mfilename('fullpath')), 'data', varargin{:});
end
