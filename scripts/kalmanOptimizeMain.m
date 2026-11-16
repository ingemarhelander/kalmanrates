% kalmanOptimizeMain.m
% Gradient ascent with backtracking line search for the Kalman objective.
% Local functions at the end of the file (requires MATLAB R2016b or later).

%% ------------------------------ Main ------------------------------------

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'line_search'));

config = read_config(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));

maxIter = config.maxIter; % Maximum number of iterations
mu = config.mu;
noisy_size = config.noisy_size; % Dimension of noisy subspace
lineSearchFcn = select_line_search(config); % Line search method and its options
dataSel = read_data_selection(config); % Days and contracts for training and test

[K, varInd, xLbd, x, data] = setup_problem(dataSel);
[f, g] = kalmanObj(x, K, mu, varInd);
history = init_history(x);
D = [];           % All orthogonal search directions 
B = [];          % Last noisy_size search directions (one per column)
evalHist = {};   % History of line search evaluations
iter = 1;
stalls = 0;      % Consecutive iterations where x did not change

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
    xStart = x;
    [x, f, g] = update_step(xNew, fNew, gNew);
    iter = iter + 1;

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
    
    %Make sure the last n search direction is stored as a base B and run noisy optimization
    [D, B] = store_direction(D, d, noisy_size);
    if iter > noisy_size
        d_prime = denoise(x, f, g, B, evalHist, noisy_size, K, mu, varInd, xLbd, config);
        if ~any(d_prime)
            fprintf('Denoising produced no feasible direction; skipping its line search.\n');
        else
            [alpha, xNew, fNew, gNew, evals] = lineSearchFcn(x, f, g, d_prime, K, mu, varInd);
            [x, f, g] = update_step(xNew, fNew, gNew);
        end

        fprintf('Denoising done \n')
    end

    % Stop if x has not changed in 3 consecutive iterations
    if isequal(x, xStart)
        stalls = stalls + 1;
        if stalls >= 3
            fprintf('Stopping: x has not changed in %d consecutive iterations.\n', stalls);
            break;
        end
    else
        stalls = 0;
    end

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
end

%% Out-of-sample test with the optimized parameters
if ~isempty(dataSel.test)
    history.outOfSample = evaluateOutOfSample(data.K, data.times, data.keepInstr, ...
        x, varInd, dataSel.train, dataSel.test);
    print_out_of_sample(history.outOfSample);
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

function [D, B] = store_direction(D, d, noisy_size)
% Keeps only the noisy_size most recent search directions (one per column).
    D = [D, d]; % adds last search direction to matrix with all previous search directions
    if size(D, 2) >= noisy_size
        B = D(:, end-noisy_size+1:end);
    else
    B = []; 
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

function [K, varInd, xLbd, x, data] = setup_problem(dataSel)
% Loads data, sets the starting point, selects which variables to optimize
% and builds the lower bounds for them. K is the training data chosen in
% the config; data holds the whole model for the out-of-sample test.

    [K, times, ricAll] = loadKalmanModel();

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

    % Parameters that belong only to contracts left out of the config do
    % not affect the likelihood, so they are not optimized (the mu penalty
    % alone would otherwise push them up without bound).
    [keepInstr, labels] = selectInstruments(ricAll, K.nu, dataSel.tenors);
    dropped = find(~keepInstr);
    droppedState = K.ns + dropped;
    unused = [thetaFInd(droppedState) SigmavInd(dropped) SigmawInd(droppedState) ...
              axInd(droppedState) Sigma0Ind(droppedState)];
    varInd = setdiff(varInd, unused, 'stable');
    if isempty(varInd)
        error('No parameters left to optimize for the chosen contracts.');
    end

    xInitAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
    xLbdAll  = [-ones(size(K.thetaF))*inf ; zeros(size(K.Sigmav)) ; ...
                 zeros(size(K.Sigmaw)) ; -inf*ones(size(K.ax)) ; ...
                 zeros(size(K.Sigma0)) ; -inf*ones(size(K.Ef(:)))];

    x    = xInitAll(varInd);
    xLbd = xLbdAll(varInd);

    % Keep the whole model (now holding the start parameters) for the test
    trainIdx = selectDays(times, dataSel.train);
    data = struct('K', K, 'times', times, 'keepInstr', keepInstr);
    K = subsetKalmanModel(K, trainIdx, keepInstr);
    print_data_selection(times(trainIdx), labels(keepInstr), dataSel);
end

function dataSel = read_data_selection(config)
% Reads the data selection keys of the config (see kalmanOptimizeMain.config).
% Missing keys mean: first/last day of the model, all weekdays, every day
% and all contracts. Without test_start there is no out-of-sample test.

    dataSel.train = struct( ...
        'startDate', config_value(config, 'train_start', []), ...
        'endDate',   config_value(config, 'train_end', []), ...
        'weekdays',  config_value(config, 'weekdays', 'all'), ...
        'dayStep',   config_value(config, 'day_step', 1));
    dataSel.tenors = config_value(config, 'tenors', 'all');
    dataSel.test = [];
    if isfield(config, 'test_start')
        dataSel.test = struct( ...
            'startDate', config.test_start, ...
            'endDate',   config_value(config, 'test_end', []));
    end
end

function value = config_value(config, key, default)
    value = default;
    if isfield(config, key)
        value = config.(key);
    end
end

function print_data_selection(trainTimes, tenorLabels, dataSel)
    fprintf('Training data: %d days, %s to %s (weekdays: %s, every %d. day)\n', ...
        numel(trainTimes), datestr(floor(trainTimes(1)), 'yyyy-mm-dd'), ...
        datestr(floor(trainTimes(end)), 'yyyy-mm-dd'), ...
        char(string(dataSel.train.weekdays)), dataSel.train.dayStep);
    fprintf('Contracts (%d): %s\n', numel(tenorLabels), strjoin(tenorLabels, ', '));
end

function print_out_of_sample(r)
    if r.conditional
        how = 'given all earlier days';
    else
        how = 'filter started at the first test day';
    end
    fprintf('\nOut-of-sample %s to %s (%s):\n', datestr(r.testStart, 'yyyy-mm-dd'), ...
        datestr(r.testEnd, 'yyyy-mm-dd'), how);
    fprintf('  Test:     log likelihood %14.4f, %4d days, per day %10.4f, per observation %8.4f\n', ...
        r.logLik, r.nDays, r.perDay, r.perObs);
    fprintf('  Training: log likelihood %14.4f, %4d days, per day %10.4f, per observation %8.4f\n', ...
        r.trainLogLik, r.trainDays, r.trainPerDay, r.trainPerObs);
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
