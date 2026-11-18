% kalmanOptimizeMain.m
% Gradient ascent for the Kalman objective, with polynomial line search
% (line_search_main) or backtracking (line_search), chosen in the config.
% Local functions at the end of the file (requires MATLAB R2016b or later).

%% ------------------------------ Main ------------------------------------

run(fullfile(fileparts(mfilename('fullpath')), '..', 'setupPaths.m'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'line_search'));

config = readConfig(fullfile(fileparts(mfilename('fullpath')), 'kalmanOptimizeMain.config'));
weights = struct('w_xi', config.w_xi, 'wbar_xi', config.wbar_xi);

maxIter = config.maxIter; % Maximum number of iterations
mu = config.mu;
noisy_size = config.noisy_size; % Dimension of noisy subspace
lineSearchFcn = select_line_search(config); % Line search method and its options
periods = experimentPeriods(config); % In-sample, training and test periods
dataSel = read_data_selection(config, periods); % Days and contracts for training and test

[K, varInd, xLbd, x, data] = setup_problem(dataSel, periods);
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
    [d, D] = find_search_direction(x, g, xLbd, D);

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
    if size(D, 2) >= noisy_size
        d_prime = denoise(x, f, g, B, evalHist, noisy_size, K, mu, varInd, xLbd, config, weights);
        if ~any(d_prime)
            fprintf('Denoising produced no feasible direction; skipping its line search.\n');
        else
            [alpha, xNew, fNew, gNew, evals] = lineSearchFcn(x, f, g, d_prime, K, mu, varInd);
            [x, f, g] = update_step(xNew, fNew, gNew);
            [D, B, evalHist] = ReOrthoginalizeAndReOrganize(d_prime, D, B, evalHist, noisy_size, evals);
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
save_result(x, varInd, data, dataSel, periods, history);

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
% Returns the line search chosen by config.lineSearch (default polynomial).
% Keys beginning with ls_ are passed, without the prefix, as options to the
% polynomial line search (see line_search_main).

    method = 'polynomial';
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
    fprintf('Main loop iteration %3d: alpha = %.4e, fNew-fPrev = %+.4e, fNew = %.6f, |gPrev| = %.4e, |gNew| = %.4e\n', ...
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


function [D, B, evalHist] = ReOrthoginalizeAndReOrganize(dprime, D, B, evalHist, noisy_size, evals)
    Dprime = ReOrthogonalize(dprime, D);
    Bprime = Dprime(:, end-noisy_size+1:end);
    evalHistPrime = ReOrganize(evalHist, evals);
    D = Dprime;
    B = Bprime;
    evalHist = evalHistPrime;
end

function Dprime = ReOrthogonalize(dprime, D)

    %Some error checks that copilot suggested
    if ~isnumeric(D) || isempty(D)
        error('D must be a nonempty numeric matrix.');
    end
    if ~isnumeric(dprime) || size(dprime, 2) ~= 1 || ...
            size(dprime, 1) ~= size(D, 1)
        error('dprime must be a numeric column vector with the same row count as D.');
    end
    if any(~isfinite(D(:))) || any(~isfinite(dprime(:)))
        error('D and dprime must contain only finite values.');
    end
    if ~isfinite(norm(dprime)) || norm(dprime) == 0
        error('dprime must have a finite, nonzero norm.');
    end
    
    %Initialize Dprime to be empty except for dprime as the last column.
    n = size(D, 2);
    Dprime = zeros(size(D), 'like', D);
    Dprime(:, end) = dprime / norm(dprime);

    % Go from n-1 down to 1 (so skip the last d that we added) and find the orthogonalized version of each column of D.
    for k = n-1:-1:1
        %extract d
        d = D(:, k);

        %Remove the components of d that are in the direction of thecolumns already in Dprime (from k+1 to n)
        for i = k+1:1:n
            d = d - (Dprime(:, i)' * d) * Dprime(:, i);
        end

        %Check that the resulting vector is finite.
        if ~isfinite(norm(d))
            error('Orthogonalization produced a non-finite residual.');
        end

        %Add in the new normalized direction.
        Dprime(:, k) = d / norm(d);
    end
end

function evalHistPrime = ReOrganize(evalHist, evals)
    %The new evalHist should be exactly evalHist except we remove the last evals that correspond to di and replace with the new evals that correspond to di prime.
    evalHistPrime = evalHist;
    evalHistPrime{end} = evals; % Rewrite the last entry with the new evaluations
end

%% Help function given by Jörgen.

function [K, varInd, xLbd, x, data] = setup_problem(dataSel, periods)
% Loads the model for the periods, sets the starting point, selects which
% variables to optimize and builds the lower bounds for them. K is the
% training data chosen in the config; data holds the whole model for the
% out-of-sample test.

    [K, times, ricAll, KI] = loadKalmanModel(periods.modelFile);
    check_model_periods(KI, periods);

    % Index ranges of each parameter block in the full parameter vector
    thetaFInd = (1:length(K.thetaF));
    SigmavInd = (1:length(K.Sigmav)) + thetaFInd(end);
    SigmawInd = (1:length(K.Sigmaw)) + SigmavInd(end);
    axInd     = (1:length(K.ax))     + SigmawInd(end);
    Sigma0Ind = (1:length(K.Sigma0)) + axInd(end);
    EfVecInd  = (1:length(K.Ef(:)))  + Sigma0Ind(end);

    % Start point: the model's own start parameters (init) or the first
    % stored optimum (saved)
    switch periods.startPoint
        case 'init'
            fprintf('Start point: the model''s own start parameters\n');
        case 'saved'
            S = load(dataPath('kalmanXoptAll.mat'), 'xOptAll');
            xOpt = S.xOptAll(:,1);
            if numel(xOpt) ~= EfVecInd(end)
                error('kalmanrates:ModelAlignment', ...
                      ['kalmanXoptAll.mat has %d parameters but the model needs %d. ' ...
                       'Use start_point init for this model.'], numel(xOpt), EfVecInd(end));
            end
            if periods.hasInsample
                warning('kalmanrates:NotOutOfSample', ...
                        ['start_point saved: kalmanXoptAll.mat was not fitted on the ' ...
                         'in-sample days only, so the test may not be out of sample.']);
            end
            K.thetaF = xOpt(thetaFInd);
            K.Sigmaw = xOpt(SigmawInd);
            K.Sigmav = xOpt(SigmavInd);
            K.ax     = xOpt(axInd);
            K.Sigma0 = xOpt(Sigma0Ind);
            K.Ef(:)  = xOpt(EfVecInd);
            fprintf('Start point: kalmanXoptAll.mat\n');
    end

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

    % Parameters that belong only to contracts left out of the config do
    % not affect the likelihood, so they are not optimized (the mu penalty
    % alone would otherwise push them up without bound).
    % Contracts without quotes on the training days cannot be fitted either.
    trainIdx = selectDays(times, dataSel.train);
    [keepInstr, labels] = selectInstruments(ricAll, K.nu, dataSel.tenors);
    quoted = false(1, K.nu);
    for t = trainIdx
        quoted = quoted | full(any(K.Iz{t}, 1));
    end
    unquoted = keepInstr & ~quoted;
    if any(unquoted)
        fprintf('Contracts without quotes on the training days are left out: %s\n', ...
                strjoin(labels(unquoted), ', '));
        keepInstr = keepInstr & quoted;
    end
    if ~any(keepInstr)
        error('kalmanrates:NoInstruments', 'None of the chosen contracts is quoted on the training days.');
    end
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
    data = struct('K', K, 'times', times, 'keepInstr', keepInstr);
    K = subsetKalmanModel(K, trainIdx, keepInstr);
    print_data_selection(times(trainIdx), labels(keepInstr), dataSel);
end

function check_model_periods(KI, periods)
% Stops when the loaded model was not built for the periods of the config.

    if ~periods.hasInsample
        return;
    end
    if ~all(isfield(KI, {'insampleStart', 'insampleEnd', 'pcaCutoff'}))
        error('kalmanrates:InvalidPeriods', ...
              '%s has no in-sample metadata. Build it again with buildKalmanModel.', periods.modelFile);
    end
    if floor(KI.insampleStart) ~= periods.insampleStart || floor(KI.insampleEnd) ~= periods.insampleEnd
        error('kalmanrates:InvalidPeriods', ...
              'The model is built on in-sample %s to %s but the config says %s to %s.', ...
              datestr(KI.insampleStart, 'yyyy-mm-dd'), datestr(KI.insampleEnd, 'yyyy-mm-dd'), ...
              datestr(periods.insampleStart, 'yyyy-mm-dd'), datestr(periods.insampleEnd, 'yyyy-mm-dd'));
    end
    if ~isempty(periods.testStart) && KI.pcaCutoff > periods.testStart
        error('kalmanrates:InvalidPeriods', ...
              'The PCA of the model uses days up to %s, after test_start %s.', ...
              datestr(KI.pcaCutoff - 1, 'yyyy-mm-dd'), datestr(periods.testStart, 'yyyy-mm-dd'));
    end
    fprintf('Model %s: in-sample %s to %s\n', periods.modelFile, ...
            datestr(KI.insampleStart, 'yyyy-mm-dd'), datestr(KI.insampleEnd, 'yyyy-mm-dd'));
end

function dataSel = read_data_selection(config, periods)
% Reads the data selection keys of the config (see kalmanOptimizeMain.config).
% The dates come from experimentPeriods. Missing keys mean: all weekdays,
% every day and all contracts. Without test_start there is no
% out-of-sample test.

    dataSel.train = struct( ...
        'startDate', periods.trainStart, ...
        'endDate',   periods.trainEnd, ...
        'weekdays',  config_value(config, 'weekdays', 'all'), ...
        'dayStep',   config_value(config, 'day_step', 1));
    dataSel.tenors = config_value(config, 'tenors', 'all');
    dataSel.test = [];
    if ~isempty(periods.testStart)
        dataSel.test = struct( ...
            'startDate', periods.testStart, ...
            'endDate',   periods.testEnd);
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

function save_result(x, varInd, data, dataSel, periods, history)
% Saves the run in data/kalmanOptimizeResult.mat for plotOutOfSample.
% xStart and xOpt are whole parameter vectors (as storeParamVec takes them).

    K = data.K;
    result.xStart = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
    result.xOpt = result.xStart;
    result.xOpt(varInd) = x;
    result.varInd = varInd;
    result.keepInstr = data.keepInstr;
    result.dataSel = dataSel;
    result.modelFile = periods.modelFile;
    result.history = history;
    result.savedAt = datestr(now, 'yyyy-mm-dd HH:MM');
    save(dataPath('kalmanOptimizeResult.mat'), '-struct', 'result');
    fprintf('Saved the result in %s\n', dataPath('kalmanOptimizeResult.mat'));
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
