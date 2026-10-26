% kalmanOptimize.m
% Gradient ascent with backtracking line search for the Kalman objective.
% Local functions at the end of the file (requires MATLAB R2016b or later).

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
    d = find_search_direction(x, g, xLbd);

    %Line search to find length of step and new X
    [alpha, xNew, fNew, gNew] = line_search(x, f, g, d, K, mu, varInd);

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
    
    %Make sure the last search direction is stored and run noisy optimization
    D = store_direction(D, d, noisy_size);
    if iter > noisy_size
        [x, f, g] = noisy_optimization(x, f, g, D, noisy_size, K, mu, varInd, xLbd);
    end

    %Check if current step fulfills stop criterion.
    if fulfills_stop_criterion(g, iter, maxIter)
        break;
    end
end

%% ------------------------------ Functions -------------------------------


function d = find_search_direction(x, g, xLbd)
% Steepest ascent direction, shortened so that a full step (alpha = 1)
% stays inside the lower bounds.

    d = g;

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


function [alpha, xNew, fNew, gNew] = line_search(x, f, g, d, K, mu, varInd)
% Backtracking line search: start with a (capped) step, halve until the
% objective improves or the maximum number of halvings is reached.

    maxLsIter = 50;

    alpha = initial_step_length(g, d);

    xNew = x + alpha*d;
    [fNew, gNew] = kalmanObj(xNew, K, mu, varInd);
    lsIter = 1;

    while true
        dfNew = gNew'*d;
        print_line_search(lsIter, alpha, fNew, f, dfNew)

        isBetter = fNew > f;
        if isBetter
            break;
        end

        alpha = alpha/2;
        xNew = x + alpha*d;
        [fNew, gNew] = kalmanObj(xNew, K, mu, varInd);

        lsIter = lsIter + 1;
        if lsIter > maxLsIter
            break;
        end
    end
end


function alpha = initial_step_length(g, d)
% Steps are far too long in the beginning, so cap the step such that the
% predicted (linear) increase in the objective is at most objInc.

    objInc   = 1E2;
    alphaMin = 1E-5;

    alpha = 1;
    alphaObjInc = objInc/(g'*d);
    alpha = min(alpha, alphaObjInc);
    alpha = max(alpha, alphaMin);
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
    fprintf('Main loop iteration %3d: alpha = %14d, fNew-fPrev = %14f, fNew = %14f, norm of gPrev = %14d, norm of gNew = %14d\n', ...
            iter, alpha, fNew-f, fNew, norm(g), norm(gNew));
end

function print_line_search(lsIter, alpha, fNew, f, dfNew)
    fprintf('LS iteration %2d: alpha = %12d, fvalNew-fval = %12f, dfNew = %12d\n', lsIter, alpha, fNew-f, dfNew);

    if isnan(fNew)
        fprintf('Warning: objective is NaN in line search\n');
    end
end

function [x, f, g] = noisy_optimization(x, f, g, D, noisy_size, K, mu, varInd, xLbd)
% Placeholder: returns x, f, g unchanged.
    fprintf('\nNoisy optimization: last %d search directions', noisy_size);
end


%% Help function we will probably not need to touch much.

%Denna har claude skrivit så kanske måste testas
function D = store_direction(D, d, noisy_size)
% Keeps only the noisy_size most recent search directions (one per column).

    D = [D, d];
    if size(D, 2) > noisy_size
        D = D(:, end-noisy_size+1:end);
    end
end

function [K, mu, maxIter, varInd, xLbd, x] = setup_problem()
% Loads data, sets the starting point, selects which variables to optimize
% and builds the lower bounds for them.

    maxIter = 100;
    % mu = 1E-4;
    % mu = 1E-8;
    mu = 1E-2;

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