% kalmanOptimizeMain.m
% Gradient ascent with backtracking line search for the Kalman objective.
% Local functions at the end of the file (requires MATLAB R2016b or later).

%% ------------------------------ Main ------------------------------------

addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'line_search'));

config = read_config(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));

maxIter = config.maxIter; % Maximum number of iterations
mu = config.mu;
noisy_size = config.noisy_size; % Dimension of noisy subspace
lineSearchFcn = select_line_search(config); % Line search method and its options

[K, varInd, xLbd, x] = setup_problem();
[f, g] = kalmanObj(x, K, mu, varInd);
history = init_history(x);
D = [];          % Last noisy_size search directions (one per column)
evalHist = {};   % History of line search evaluations
iter = 1;

while true
    history = update_history(history, x, f, g);

    %Find search direction d
    d = find_search_direction(x, g, xLbd, D);

    %Line search to find length of step and new X
    [alpha, xNew, fNew, gNew, evals] = lineSearchFcn(x, f, g, d, K, mu, varInd);

    %Print and update history
    print_iteration(iter, alpha, f, fNew, g, gNew);
    history.stepSizes(end+1, 1) = alpha;
    evalHist = store_evals(evalHist, evals, noisy_size); 

    %Iteration done, update counter and x f and g
    [x, f, g] = update_step(xNew, fNew, gNew);
    iter = iter + 1;

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
    
    %Make sure the last search direction is stored and run noisy optimization
    D = store_direction(D, d, noisy_size);
    if iter > noisy_size
        [x, f, g] = denoise(x, f, g, D, evalHist, noisy_size, K, mu, varInd, xLbd);
    end

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
end

%% ------------------------------ Functions -------------------------------

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

function lineSearchFcn = select_line_search(config)
% Returns the line search chosen by config.lineSearch (default backtracking).
% Keys beginning with ls_ are passed, without the prefix, as options to the
% polynomial line search (see line_search_main).

    method = 'backtracking';
    if isfield(config, 'lineSearch')
        method = config.lineSearch;
    end

    lsOpts = struct();
    keys = fieldnames(config);
    for i = 1:numel(keys)
        if startsWith(keys{i}, 'ls_')
            lsOpts.(keys{i}(4:end)) = config.(keys{i});
        end
    end

    switch lower(method)
        case 'polynomial'
            lineSearchFcn = @(x, f, g, d, K, mu, varInd) ...
                line_search_main(x, f, g, d, K, mu, varInd, lsOpts);
        case 'backtracking'
            if ~isempty(fieldnames(lsOpts))
                warning('ls_ options are ignored by the backtracking line search.');
            end
            lineSearchFcn = @line_search;
        otherwise
            error('Unknown lineSearch in config: %s (use polynomial or backtracking)', method);
    end
    fprintf('Line search: %s\n', lower(method));
end

function print_iteration(iter, alpha, f, fNew, g, gNew)
    fprintf('Main loop iteration %3d: alpha = %14d, fNew-fPrev = %14f, fNew = %14f, norm of gPrev = %14d, norm of gNew = %14d\n', ...
            iter, alpha, fNew-f, fNew, norm(g), norm(gNew));
end

function D = store_direction(D, d, noisy_size)
% Keeps only the noisy_size most recent search directions (one per column).

    D = [D, d];
    if size(D, 2) > noisy_size
        D = D(:, end-noisy_size+1:end);
    end
end

function evalHist = store_evals(evalHist, evals, noisy_size)
% Keeps the evaluations (x, f, g) of the noisy_size most recent steps.
% evalHist{k} belongs to step k of the window (oldest first), so it lines
% up with the columns of D.

    evalHist{end+1} = evals;
    if numel(evalHist) > noisy_size
        evalHist = evalHist(end-noisy_size+1:end);
    end
end

%% Help function given by Jörgen.

function [K, varInd, xLbd, x] = setup_problem()
% Loads data, sets the starting point, selects which variables to optimize
% and builds the lower bounds for them.

    S = load('kalmanDebug.mat');
    if isfield(S, 'K')
        K = S.K;
    else
        S = load('kalmanDebugR.mat');
        K = S.K;
    end

    % Index ranges of each parameter block in the full parameter vector
    thetaFInd = (1:length(K.thetaF));
    SigmavInd = (1:length(K.Sigmav)) + thetaFInd(end);
    SigmawInd = (1:length(K.Sigmaw)) + SigmavInd(end);
    axInd     = (1:length(K.ax))     + SigmawInd(end);
    Sigma0Ind = (1:length(K.Sigma0)) + axInd(end);
    EfVecInd  = (1:length(K.Ef(:)))  + Sigma0Ind(end);

    % Start from the first stored optimum
    S = load('kalmanXoptAll.mat');
    xOpt = S.xOptAll(:,1);
    K.thetaF = xOpt(thetaFInd);
    K.Sigmaw = xOpt(SigmawInd);
    K.Sigmav = xOpt(SigmavInd);
    K.ax     = xOpt(axInd);
    K.Sigma0 = xOpt(Sigma0Ind);
    K.Ef(:)  = xOpt(EfVecInd);

    % Which variables to optimize
    optCase = 8;
    switch optCase
        case 1, varInd = [thetaFInd SigmavInd SigmawInd axInd Sigma0Ind EfVecInd];
        case 2, varInd = [SigmavInd SigmawInd axInd Sigma0Ind];
        case 3, varInd = [SigmavInd SigmawInd Sigma0Ind];
        case 4, varInd = [EfVecInd];
        case 5, varInd = [SigmavInd];
        case 6, varInd = [thetaFInd SigmavInd SigmawInd Sigma0Ind];
        case 7, varInd = [thetaFInd];
        case 8, varInd = [SigmawInd];
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

function config = read_config(configPath)
% Reads key/value rows into a struct. Blank lines and lines beginning with
% # or % are ignored.

    fid = fopen(configPath, 'r');
    if fid == -1
        error('Could not open config file: %s', configPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    config = struct();

    while true
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        line = strtrim(line);
        if isempty(line) || startsWith(line, '#') || startsWith(line, '%')
            continue;
        end

        tokens = regexp(line, '^([A-Za-z]\w*)\s+(.+)$', 'tokens', 'once');
        if isempty(tokens) || ~isvarname(tokens{1})
            error('Invalid config row: %s', line);
        end

        key = tokens{1};
        valueText = strtrim(tokens{2});
        if strcmpi(valueText, 'true')
            value = true;
        elseif strcmpi(valueText, 'false')
            value = false;
        else
            numericValue = str2double(valueText);
            if ~isnan(numericValue) || strcmpi(valueText, 'nan')
                value = numericValue;
            else
                value = valueText;
            end
        end
        config.(key) = value;
    end
end