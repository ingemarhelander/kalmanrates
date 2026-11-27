function kalmanOptimizeSDLineSearch(config)
% Gradient ascent. Defaults to backtracking; configure lineSearch as
% 'polynomial' or 'backtracking', with optional ls_* polynomial options.
% Local functions at the end of the file (requires MATLAB R2016b or later).

    if nargin < 1 || isempty(config)
        config = struct();
    end
    if ischar(config) || (isstring(config) && isscalar(config))
        config = struct('lineSearch', char(config));
    end
    if ~isstruct(config) || ~isscalar(config)
        error('kalmanOptimizeSDLineSearch:InvalidConfig', ...
              'Input must be a scalar config struct or a line-search method name.');
    end
    if ~isfield(config, 'lineSearch')
        config.lineSearch = 'backtracking';
    end
    if isstring(config.lineSearch) && isscalar(config.lineSearch)
        config.lineSearch = char(config.lineSearch);
    end
    if ~ischar(config.lineSearch) || size(config.lineSearch, 1) ~= 1
        error('kalmanOptimizeSDLineSearch:InvalidLineSearch', ...
              'config.lineSearch must be ''backtracking'' or ''polynomial''.');
    end

    run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
    lineSearchFcn = select_line_search(config);

%% ------------------------------ Main ------------------------------------
[K, mu, maxIter, varInd, xLbd, x] = setup_problem();
noisy_size = 3; % Dimension of noisy subspace
[f, g] = kalmanObj(x, K, mu, varInd);
history = init_history(x);
D = [];          % Last noisy_size search directions (one per column)
iter = 1;

while true
    history = update_history(history, x, f, g);

    %Find search direction d
    d = find_search_direction(x, g, xLbd, strcmpi(config.lineSearch, 'polynomial'));

    %Line search to find length of step and new X
    [alpha, xNew, fNew, gNew] = lineSearchFcn(x, f, g, d, K, mu, varInd);

    %Print and update history
    print_iteration(iter, alpha, f, fNew, g, gNew);
    history.stepSizes(end+1, 1) = alpha;

    %Iteration done, update counter and x f and g
    [x, f, g] = update_step(xNew, fNew, gNew);
    iter = iter + 1;

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
  
    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
end
end

%% ------------------------------ Functions -------------------------------

function lineSearchFcn = select_line_search(config)
% Selects the shared line search and passes ls_* fields as polynomial options.

    lineSearchPath = fullfile(fileparts(mfilename('fullpath')), 'line_search');
    addpath(lineSearchPath);

    method = lower(config.lineSearch);
    lsOpts = struct();
    keys = fieldnames(config);
    for i = 1:numel(keys)
        if startsWith(keys{i}, 'ls_')
            lsOpts.(keys{i}(4:end)) = config.(keys{i});
        end
    end

    switch method
        case 'polynomial'
            if ~isfield(lsOpts, 'alphaMin')
                % This script's gradient scale needs smaller initial steps.
                lsOpts.alphaMin = 1e-16;
            end
            lineSearchFcn = @(x, f, g, d, K, mu, varInd) ...
                line_search_main(x, f, g, d, K, mu, varInd, lsOpts);
        case 'backtracking'
            if ~isempty(fieldnames(lsOpts))
                warning('kalmanrates:IgnoredLineSearchOptions', ...
                        'ls_ options are ignored by the backtracking line search.');
            end
            lineSearchFcn = @line_search;
        otherwise
            error('kalmanOptimizeSDLineSearch:InvalidLineSearch', ...
                  'Unknown lineSearch in config: %s (use polynomial or backtracking).', ...
                  config.lineSearch);
    end
    fprintf('Line search: %s\n', method);
end


function d = find_search_direction(x, g, xLbd, normalizeDirection)
% Steepest ascent direction, shortened so that a full step (alpha = 1)
% stays inside the lower bounds.

    d = g;
    if normalizeDirection && norm(d) > 0
        d = d / norm(d);
    end

    if g'*d < 0
        fprintf('This is not an ascent direction\n');
    end

    % Largest step length (along d) before each variable hits its lower bound
    alphaLbd = ones(size(x))*inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);

    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d*alphaLbdMin*0.99;   % Ensures that alpha = 1 is feasible
    end

    if any(isnan(d))
        error('Search direction contains NaN');
    end
end

function stop = fulfills_stop_criterion(g, iter, maxIter)
% Stop when the gradient is small or the iteration limit is reached.

    stop = false;
    if norm(g) < 1e-7
        fprintf('The norm of the gradient: %+.5e\n', norm(g));
        stop = true;
    elseif iter >= maxIter
        fprintf('Maximum number of iterations\n');
        stop = true;
    end
end


function print_iteration(iter, alpha, f, fNew, g, gNew)
    fprintf('Main loop iteration %3d: alpha = %.3e, fNew-fPrev = %+.6e, fNew = %.6f, norm of gPrev = %.3e, norm of gNew = %.3e\n', ...
            iter, alpha, fNew-f, fNew, norm(g), norm(gNew));
end

%% Help function we will probably not need to touch much.

function [K, mu, maxIter, varInd, xLbd, x] = setup_problem()
% Loads data, sets the starting point, selects which variables to optimize
% and builds the lower bounds for them.

    maxIter = 100;
    % mu = 1E-4;
    % mu = 1E-8;
    mu = 1E-2;

    [K, ~] = loadKalmanSample();

    % Index ranges of each parameter block in the full parameter vector
    thetaFInd = (1:length(K.thetaF));
    SigmavInd = (1:length(K.Sigmav)) + thetaFInd(end);
    SigmawInd = (1:length(K.Sigmaw)) + SigmavInd(end);
    axInd     = (1:length(K.ax))     + SigmawInd(end);
    Sigma0Ind = (1:length(K.Sigma0)) + axInd(end);
    EfVecInd  = (1:length(K.Ef(:)))  + Sigma0Ind(end);

    % Start from the first stored optimum
    S = load(dataPath('kalmanXoptAll.mat'), 'xOptAll');
    xOpt = S.xOptAll(:,1);
    K.thetaF = xOpt(thetaFInd);
    K.Sigmaw = xOpt(SigmawInd);
    K.Sigmav = xOpt(SigmavInd);
    K.ax     = xOpt(axInd);
    K.Sigma0 = xOpt(Sigma0Ind);
    K.Ef(:)  = xOpt(EfVecInd);

    % Which variables to optimize
    optCase = 9;
    switch optCase
        case 1, varInd = [thetaFInd SigmavInd SigmawInd axInd Sigma0Ind EfVecInd];
        case 2, varInd = [SigmavInd SigmawInd axInd Sigma0Ind];
        case 3, varInd = [SigmavInd SigmawInd Sigma0Ind];
        case 4, varInd = [EfVecInd];
        case 5, varInd = [SigmavInd];
        case 6, varInd = [thetaFInd SigmavInd SigmawInd Sigma0Ind];
        case 7, varInd = [thetaFInd];
        case 8, varInd = [SigmawInd];
        case 9, varInd = [thetaFInd SigmavInd SigmawInd axInd Sigma0Ind];
    end

    xInitAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
    xLbdAll  = [-ones(size(K.thetaF))*inf ; zeros(size(K.Sigmav)) ; ...
                 zeros(size(K.Sigmaw)) ; -inf*ones(size(K.ax)) ; ...
                 zeros(size(K.Sigma0)) ; -inf*ones(size(K.Ef(:)))];

    x    = xInitAll(varInd);
    xLbd = xLbdAll(varInd);
end

function [x, f, g] = update_step(xNew, fNew, gNew)
% Take the step.

    x = xNew;
    f = fNew;
    g = gNew;
end

function history = init_history(xInit)
    history.xInit = xInit;
    history.fvals = [];
    history.normOfgrads = [];
    history.normOftravelledDistances = [];
    history.stepSizes = [];
    history.xs = [];
    history.grads = [];
end


function history = update_history(history, x, f, g)
% Records the state at the start of each iteration.

    history.fvals(end+1, 1) = f;
    history.normOfgrads(end+1, 1) = norm(g);
    history.normOftravelledDistances(end+1, 1) = norm(history.xInit - x);
    history.xs(:, end+1) = x;
    history.grads(:, end+1) = g;
end