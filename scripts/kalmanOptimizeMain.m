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
periods = experimentPeriods(config); % In-sample, optimization and out-of-sample periods

[K, varInd, xLbd, x, data] = setup_problem(periods);
lineSearchFcn = select_line_search(config, xLbd); % Line search method and its options
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

%% Out-of-sample evaluation with the optimized parameters
if ~isempty(periods.outOfSampleStart)
    history.outOfSample = evaluateOutOfSample(data.K, data.times, x, varInd, periods);
    print_out_of_sample(history.outOfSample);
end
save_result(x, varInd, data.K, periods, history);

% Full parameter vector for runCurvePipelinePlots. K is a struct here (from
% subsetKalmanModel), so the vector is built from the parameter fields.
parameterVector = [K.thetaF; K.Sigmav; K.Sigmaw; K.ax; K.Sigma0; K.Ef(:)];
parameterVector(varInd) = x;
kalmanOptimizeMainOutput.parameterVector = parameterVector;
save(dataPath('kalmanOptimizeMainOutput.mat'), 'kalmanOptimizeMainOutput');

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

function lineSearchFcn = select_line_search(config, xLbd)
% Returns the line search chosen by config.lineSearch (default polynomial).
% Keys beginning with ls_ are passed, without the prefix, as options to the
% polynomial line search (see line_search_main). The lower bounds xLbd are
% passed too, so that no step goes past the first variance bound along d.

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
            lsOpts.xLbd = xLbd;
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

function [K, varInd, xLbd, x, data] = setup_problem(periods)
% Loads the model for the periods, sets the starting point, selects which
% variables to optimize and builds the lower bounds for them. K holds the
% optimization days; data holds the whole model for the out-of-sample
% evaluation.

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
                         'in-sample days only, so the evaluation may not be out of sample.']);
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

    % Contracts without quotes on the optimization days do not affect the
    % likelihood, so their parameters are not optimized (the mu penalty
    % alone would otherwise push them up without bound).
    optimizeIdx = selectDays(times, periods.optimizeStart, periods.optimizeEnd);
    quoted = false(1, K.nu);
    for t = optimizeIdx
        quoted = quoted | full(any(K.Iz{t}, 1));
    end
    if ~all(quoted)
        unquoted = find(~quoted);
        fprintf('Contracts without quotes on the optimization days are not optimized: %s\n', ...
                strjoin(ricAll(unquoted), ', '));
        state = K.ns + unquoted;
        varInd = setdiff(varInd, [thetaFInd(state) SigmavInd(unquoted) SigmawInd(state) ...
                                  axInd(state) Sigma0Ind(state)], 'stable');
    end

    xInitAll = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
    xLbdAll  = [-ones(size(K.thetaF))*inf ; zeros(size(K.Sigmav)) ; ...
                 zeros(size(K.Sigmaw)) ; -inf*ones(size(K.ax)) ; ...
                 zeros(size(K.Sigma0)) ; -inf*ones(size(K.Ef(:)))];

    x    = xInitAll(varInd);
    xLbd = xLbdAll(varInd);

    % Keep the whole model (now holding the start parameters) for the
    % out-of-sample evaluation
    data = struct('K', K, 'times', times);
    K = subsetKalmanModel(K, optimizeIdx);
    fprintf('Optimization data: %d days, %s to %s\n', numel(optimizeIdx), ...
            datestr(floor(times(optimizeIdx(1))), 'yyyy-mm-dd'), datestr(floor(times(optimizeIdx(end))), 'yyyy-mm-dd'));
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
    if ~isempty(periods.outOfSampleStart) && KI.pcaCutoff > periods.outOfSampleStart
        error('kalmanrates:InvalidPeriods', ...
              'The PCA of the model uses days up to %s, after out_of_sample_start %s.', ...
              datestr(KI.pcaCutoff - 1, 'yyyy-mm-dd'), datestr(periods.outOfSampleStart, 'yyyy-mm-dd'));
    end
    fprintf('Model %s: in-sample %s to %s\n', periods.modelFile, ...
            datestr(KI.insampleStart, 'yyyy-mm-dd'), datestr(KI.insampleEnd, 'yyyy-mm-dd'));
end

function print_out_of_sample(r)
    fprintf('\nOut-of-sample %s to %s (given all earlier days):\n', ...
        datestr(r.outOfSampleStart, 'yyyy-mm-dd'), datestr(r.outOfSampleEnd, 'yyyy-mm-dd'));
    fprintf('  Out-of-sample: log likelihood %14.4f, %4d days, per day %10.4f, per observation %8.4f\n', ...
        r.logLik, r.nDays, r.perDay, r.perObs);
    fprintf('  Optimization:  log likelihood %14.4f, %4d days, per day %10.4f, per observation %8.4f\n', ...
        r.optimizeLogLik, r.optimizeDays, r.optimizePerDay, r.optimizePerObs);
end

function save_result(x, varInd, K, periods, history)
% Saves the run in data/kalmanOptimizeResult.mat for plotOutOfSample.
% K is the whole model holding the start parameters. xStart and xOpt are
% whole parameter vectors (as storeParamVec takes them).

    result.xStart = [K.thetaF ; K.Sigmav ; K.Sigmaw ; K.ax ; K.Sigma0 ; K.Ef(:)];
    result.xOpt = result.xStart;
    result.xOpt(varInd) = x;
    result.varInd = varInd;
    result.periods = periods;
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
