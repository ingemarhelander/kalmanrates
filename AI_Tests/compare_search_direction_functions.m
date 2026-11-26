function compare_search_direction_functions()
% Compare the original and orthogonal search-direction implementations.

    root = fileparts(fileparts(mfilename('fullpath')));
    functionPath = fullfile(root, 'AI_Implementations');
    originalFunctionPath = fullfile(root, 'scripts', 'functions');
    addpath(functionPath);
    addpath(originalFunctionPath);
    pathCleanup = onCleanup(@() rmpath(functionPath)); %#ok<NASGU>
    originalPathCleanup = onCleanup(@() rmpath(originalFunctionPath)); %#ok<NASGU>

    cases = {
        'Initial direction', ...
        [0; 0], [3; 4], [-inf; -inf], [];
        'Existing orthogonal history', ...
        zeros(3, 1), [2; 3; 4], -inf(3, 1), [1; 0; 0];
        'Lower-bound handling', ...
        [0.1; 2], [-2; 1], [0; -inf], [];
        'Gradient in history span', ...
        zeros(3, 1), [3; 4; 0], -inf(3, 1), [0.6; 0.8; 0]
    };

    fprintf('%-30s %-12s %-12s %-14s %-14s\n', ...
        'Case', 'Old norm', 'New norm', 'Directions equal', ...
        'Histories equal');
    fprintf('%s\n', repmat('-', 1, 84));

    for k = 1:size(cases, 1)
        name = cases{k, 1};
        x = cases{k, 2};
        g = cases{k, 3};
        xLbd = cases{k, 4};
        D = cases{k, 5};

        [dOld, Dold] = find_search_direction(x, g, xLbd, D);
        [dNew, Dnew] = find_search_direction_orthogonal(x, g, xLbd, D);

        sameDirection = norm(dOld - dNew) <= 1e-12;
        sameHistory = isequal(Dold, Dnew);

        fprintf('%-30s %-12.4g %-12.4g %-14s %-14s\n', ...
            name, norm(dOld), norm(dNew), ...
            string(sameDirection), string(sameHistory));

        assert(all(isfinite(dOld)) && all(isfinite(dNew)), ...
            'Both functions must return finite directions.');
        assert(all(x + dOld >= xLbd) && all(x + dNew >= xLbd), ...
            'Both directions must respect lower bounds for a unit step.');
        assert(sameDirection && sameHistory, ...
            'Both implementations should return the same direction and history.');
        if ~isempty(Dnew)
            assert(norm(Dnew' * dNew) <= 1e-12, ...
                'The new direction must be orthogonal to retained history.');
        end
    end

    fprintf('\nAll comparison checks passed.\n');
end
