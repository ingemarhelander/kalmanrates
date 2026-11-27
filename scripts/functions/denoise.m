function d_prime = denoise(x, f, g, B, evalHist, noisy_size, K, mu, varInd, xLbd, config, weights)
% Placeholder: returns x, f, g unchanged.
    B = orth(B);  % Ensure B is orthonormal
    columnNorms = sqrt(sum(abs(B).^2, 1));

    %Check that no columns are of length 0.
    if any(columnNorms == 0)
        error('denoise:ZeroColumnInD', 'B cannot contain zero columns.');
    end

    % Normalize the columns of B to get an orthonormal basis for the span of the past search directions.
    B = B ./ columnNorms;

    % Get the X, F and G values from evalHist.
    X = [];
    F = [];
    G = [];
    for k = 1:numel(evalHist)
        X = [X, evalHist{k}.X];
        F = [F, evalHist{k}.F];
        G = [G, evalHist{k}.G];
    end
    % Keep the current evaluation last so the model is centered at x.
    currentMatches = all(X == x, 1) & (F == f) & all(G == g, 1);
    X(:, currentMatches) = [];
    F(currentMatches) = [];
    G(:, currentMatches) = [];
    X = [X, x];
    F = [F, f];
    G = [G, g];
    % Eval j corresponds to X(:,j), F(j) and G(:,j)

    checkSpread = true;
    if isfield(config, 'denoise_check_spread')
        checkSpread = config.denoise_check_spread;
    end

    spreadRemovalMethod = 0;
    if isfield(config, 'denoise_spread_removal_method')
        spreadRemovalMethod = config.denoise_spread_removal_method;
    end

    distanceWeightScale = 1;
    if isfield(config, 'denoise_distance_weight_scale')
        distanceWeightScale = config.denoise_distance_weight_scale;
    end

    [Q, q, c, e, ok] = get_quadratic_approx_in_subspace( ...
        X, F, G, B, config.fixed_or_exact, weights, ...
        logical(checkSpread), double(spreadRemovalMethod), distanceWeightScale);
    if ~ok
        d_prime = zeros(size(x));
        return;
    end

    % Check that q is finite.
    if any(~isfinite(q))
        error('denoise:NonFiniteDirection', 'The quadratic model produced a non-finite direction.');
    end

    % Compute the search direction d_prime from the quadratic model.
    qNorm = norm(q);
    if qNorm == 0
        d_prime = zeros(size(q));
        return;
    end
    d_prime = q / qNorm;

    % Ascent direction: if the cosine of the angle between d_prime and g is
    % below minCos, rotate d_prime towards g in the plane spanned by d_prime
    % and g, just far enough that the cosine is exactly minCos. The result
    % is the first such direction on the convex combinations from d_prime
    % to g/||g||, normalized.
    minCos = 1e-7;
    if isfield(config, 'denoise_min_cos_angle')
        minCos = config.denoise_min_cos_angle;
    end

    %if g' * d_prime <= 0 && any(g), d_prime = 0.5 * (d_prime + g / norm(g)); end


    %{}

    if any(g)
        gHat = g / norm(g);
        if gHat' * d_prime < minCos
            p = d_prime - (gHat' * d_prime) * gHat;   % part of d_prime orthogonal to g
            if norm(p) > 0
                d_prime = minCos * gHat + sqrt(1 - minCos^2) * p / norm(p);
            else
                d_prime = gHat;                       % d_prime = -gHat: no plane, use g
            end
        end
    end
    %}

    % Make sure tyhe search direction d_prime is feasible with respect to the lower bounds xLbd.
    alphaLbd = ones(size(x)) * inf;
    hitsBound = (d_prime < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d_prime(hitsBound);
    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d_prime = d_prime * alphaLbdMin * 0.99;
    end
end


function [Q, q, c, e, ok] = get_quadratic_approx_in_subspace( ...
    X, F, G, B, method, weights, checkSpread, spreadRemovalMethod, ...
    distanceWeightScale)
% Build S_bar for the full basis first. Only inspect sample spread and retry
% in a smaller subspace when the full-basis system cannot be used reliably.
    originalDimension = size(B, 1);
    originalB = B;
    initialBasisDimension = size(B, 2);
    removedDimensions = 0;
    Z = B' * X;
    Gsmall = B' * G;

    [modelData, ok, illConditioned, failureReason, regularization] = ...
        prepare_quadratic_approximation(Z, weights, distanceWeightScale);
    originalModelData = modelData;
    originalOk = ok;
    originalIllConditioned = illConditioned;
    originalFailureReason = failureReason;
    originalRegularization = regularization;

    if (~ok || illConditioned) && checkSpread
        [reducedB, removedCount] = ...
            remove_poorly_spread_directions(X, B, spreadRemovalMethod);
        if removedCount > 0
            removedDimensions = removedCount;
            if size(reducedB, 2) == 0
                ok = false;
                failureReason = 'spread check removed all basis directions';
                if originalOk
                    B = originalB;
                    modelData = originalModelData;
                    ok = true;
                    illConditioned = originalIllConditioned;
                    removedDimensions = 0;
                    failureReason = '';
                end
            else
                reducedZ = reducedB' * X;
                reducedGsmall = reducedB' * G;
                [reducedModelData, reducedOk, reducedIllConditioned, ...
                    reducedFailureReason, reducedRegularization] = ...
                    prepare_quadratic_approximation( ...
                        reducedZ, weights, distanceWeightScale);
                if reducedOk
                    B = reducedB;
                    Z = reducedZ;
                    Gsmall = reducedGsmall;
                    modelData = reducedModelData;
                    ok = true;
                    illConditioned = reducedIllConditioned;
                    failureReason = '';
                    regularization = reducedRegularization;
                elseif originalOk
                    B = originalB;
                    modelData = originalModelData;
                    ok = true;
                    illConditioned = originalIllConditioned;
                    removedDimensions = 0;
                    failureReason = '';
                    regularization = originalRegularization;
                else
                    B = reducedB;
                    Z = reducedZ;
                    Gsmall = reducedGsmall;
                    modelData = reducedModelData;
                    ok = false;
                    illConditioned = false;
                    failureReason = reducedFailureReason;
                    regularization = 0;
                end
            end
        elseif ~ok
            failureReason = sprintf('%s; spread check found no removable directions', ...
                originalFailureReason);
        end
    end

    removedDimensions = initialBasisDimension - size(B, 2);
    if regularization > 0
        regularizationStatus = sprintf('yes (s=%.3e)', regularization);
    else
        regularizationStatus = 'no';
    end
    if ok
        if illConditioned
            conditionStatus = 'ill-conditioned';
        else
            conditionStatus = 'well-conditioned';
        end
        fprintf(['[denoise] Using %d of %d dimensions (%d removed); ' ...
            'S_bar is %s; regularization applied: %s.\n'], ...
            size(B, 2), initialBasisDimension, removedDimensions, ...
            conditionStatus, regularizationStatus);
    else
        fprintf(['[denoise] No model computed; 0 of %d dimensions used ' ...
            '(%d removed); conditioning unavailable. Regularization applied: %s. ' ...
            'Reason: %s.\n'], initialBasisDimension, removedDimensions, ...
            regularizationStatus, failureReason);
    end

    if ~ok
        Q = NaN(originalDimension, originalDimension);
        q = NaN(originalDimension, 1);
        c = NaN;
        e = NaN(originalDimension + 1, size(X, 2));
        return;
    end

    [Qsmall, qsmall, c, esmall] = ...
        compute_quadratic_approximation(F, Gsmall, modelData, method);

    Q = B * Qsmall * B';
    q = B * qsmall;
    e = [esmall(1, :); B * esmall(2:end, :)];
    ok = true;
end


function [B, removedCount] = remove_poorly_spread_directions(X, B, spreadRemovalMethod)
% Reduce B using the configured spread criterion after full S_bar fails.
    if ~ismember(spreadRemovalMethod, [0, 1])
        error('denoise:InvalidRemovingOption', ...
            'Spread removal method must be 0 (B columns) or 1 (SVD directions).');
    end

    centeredCoordinates = B' * X;
    centeredCoordinates = centeredCoordinates - mean(centeredCoordinates, 2);
    if spreadRemovalMethod == 0
        directionSpread = sqrt(sum(centeredCoordinates.^2, 2));
        directionVectors = eye(size(B, 2));
    else
        [spreadDirections, singularValues, ~] = ...
            svd(centeredCoordinates, 'econ');
        directionSpread = diag(singularValues);
        directionSpread(end + 1:size(B, 2), 1) = 0;
        directionVectors = spreadDirections;
    end

    minimumSpreadRatio = sqrt(1e-12);
    if isempty(directionSpread)
        maximumSpread = 0;
    else
        maximumSpread = max(directionSpread);
    end
    minimumSpread = minimumSpreadRatio * maximumSpread;
    if maximumSpread == 0
        retained = false(size(directionSpread));
    else
        retained = directionSpread > minimumSpread;
    end

    removedCount = nnz(~retained);
    if any(retained)
        if spreadRemovalMethod == 0
            B = B(:, retained);
        else
            B = B * directionVectors(:, retained(1:size(directionVectors, 2)));
        end
    else
        B = B(:, false(1, size(B, 2)));
    end
end


function [model, ok, illConditioned, failureReason, regularization] = ...
    prepare_quadratic_approximation(X, weights, distanceWeightScale)
% Build, scale, and condition-check S_bar and the data used by the fit.
    model = struct();
    ok = false;
    illConditioned = false;
    failureReason = '';
    regularization = 0;
    n = size(X, 1);
    m = size(X, 2);
    deltaX = X - X(:, end);
    nVech = n * (n + 1) / 2;
    p = nVech + n + 1;

    distances = sqrt(sum(deltaX.^2, 1));
    meanDistance = mean(distances);
    if meanDistance == 0
        distanceWeights = ones(1, m);
    else
        scaledDistances = distances ./ meanDistance ./ distanceWeightScale;
        distanceWeights = exp(-0.5 * scaledDistances.^2);
    end

    w_xi = (ones(1, m) * weights.w_xi) .* distanceWeights;
    wbar_xi = (ones(1, m) * weights.wbar_xi) .* distanceWeights;
    w = 1 + w_xi;
    wbar = 1 + wbar_xi;

    duplication = zeros(n * n, nVech);
    k = 0;
    for col = 1:n
        for row = col:n
            k = k + 1;
            idxLower = (col - 1) * n + row;
            idxUpper = (row - 1) * n + col;
            duplication(idxLower, k) = 1;
            duplication(idxUpper, k) = 1;
        end
    end

    a_all = zeros(p, m);
    A_all = zeros(n, p, m);
    for i = 1:m
        dx = deltaX(:, i);
        halfOuter = 0.5 * (dx * dx');
        topPart = duplication' * halfOuter(:);
        a_all(:, i) = [topPart; dx; 1];
        A_all(:, :, i) = [kron(dx', eye(n)) * duplication, eye(n), zeros(n, 1)];
    end

    S_bar = zeros(p, p);
    for i = 1:m
        a_i = a_all(:, i);
        A_i = A_all(:, :, i);
        S_bar = S_bar + w_xi(i) * (a_i * a_i') ...
            + wbar_xi(i) * (A_i' * A_i);
    end
    S_bar = 0.5 * (S_bar + S_bar');
    if any(~isfinite(S_bar(:)))
        failureReason = 'S_bar contains nonfinite entries';
        return;
    end
    diagonal = diag(S_bar);
    scale = ones(p, 1);
    positiveDiagonal = diagonal > 0;
    scale(positiveDiagonal) = 1 ./ sqrt(diagonal(positiveDiagonal));
    a_all = scale .* a_all;
    A_all = A_all .* scale';
    S_bar = S_bar .* (scale * scale');
    [R, notPosDef] = chol(S_bar);
    if notPosDef
        epsilon = 1e-12; % Minimum eigenvalue floor for the scaled matrix.
        lambdaMin = min(eig(S_bar));
        regularization = max(0, epsilon - lambdaMin);
        S_bar = S_bar + regularization * eye(p);
        [R, notPosDef] = chol(S_bar);
    end
    if notPosDef
        failureReason = 'Cholesky factorization of scaled S_bar failed';
        return;
    end

    reciprocalCondition = rcond(S_bar);
    illConditioned = reciprocalCondition < 1e-12;

    model = struct('n', n, 'm', m, 'nVech', nVech, 'p', p, ...
        'w_xi', w_xi, 'wbar_xi', wbar_xi, 'w', w, 'wbar', wbar, ...
        'a_all', a_all, 'A_all', A_all, 'scale', scale, ...
        'duplication', duplication, 'R', R);
    ok = true;
end


function [Q, q, c, e] = compute_quadratic_approximation(F, G, model, method)
% Fit the quadratic model using prepared data and the checked S_bar factor.

    n = model.n;
    m = model.m;
    nVech = model.nVech;
    p = model.p;
    w_xi = model.w_xi;
    wbar_xi = model.wbar_xi;
    w = model.w;
    wbar = model.wbar;
    a_all = model.a_all;
    A_all = model.A_all;
    s = model.scale;
    duplication = model.duplication;
    R = model.R;

    % Factorise once. S_bar is constant, so the factor can be reused in
    % every iteration. Needs S_bar positive definite, which requires
    % enough points in general position (see note below).
    %   S_bar = R' * R,  so  S_bar \ b = R \ (R' \ b)

    % ------------------------------------------------------------------
    % Steps 5-7: solve the coupled equations (96)-(97) for the noise
    % ------------------------------------------------------------------
    % Returns e_f (1 x m) and e_g (n x m).

    if method == 0
        [e_f, e_g] = fixed_point(F, G, ...
                                 a_all, A_all, ...
                                 w_xi, wbar_xi, w, wbar, ...
                                 R);
    elseif method == 1
        [e_f, e_g] = solve_noise_system(F, G, ...
                                        a_all, A_all, ...
                                        w_xi, wbar_xi, w, wbar, ...
                                        R);
    else
        error('denoise:UnknownMethod', 'Unknown method for solving the noise system.');
    end

    % ------------------------------------------------------------------
    % Recompute y from the final noise terms (Step 3 and Step 4)
    % ------------------------------------------------------------------
    % Noise-corrected data: r_i = f_i - e_i,  rbar_i = g_i - ebar_i
    r    = F - e_f;            % 1 x m
    rbar = G - e_g;            % n x m

    % alpha = sum_i ( w_xi_i * r_i * a_i + wbar_xi_i * A_i' * rbar_i )   (S3b)
    alpha = zeros(p, 1);

    for i = 1:m
        a_i = a_all(:, i);
        A_i = A_all(:, :, i);

        alpha_function = w_xi(i) * r(i) * a_i;
        alpha_gradient = wbar_xi(i) * (A_i' * rbar(:, i));

        alpha = alpha + alpha_function + alpha_gradient;
    end

    % Best y for these noise terms:  y = S_bar^{-1} alpha   (86)
    y = R \ (R' \ alpha);
    y = s .* y;

    % ------------------------------------------------------------------
    % Unpack y = [vech(Q); q; c] into Q, q, c
    % ------------------------------------------------------------------
    vechQ = y(1 : nVech);
    q     = y(nVech + 1 : nVech + n);
    c     = y(end);

    % vec(Q) = duplication * vech(Q), then reshape into an n x n matrix.
    vecQ = duplication * vechQ;
    Q    = reshape(vecQ, n, n);

    % Remove tiny asymmetry from rounding (Q is symmetric in theory)
    Q = 0.5 * (Q + Q');

    % ------------------------------------------------------------------
    % Pack the noise terms into one output matrix
    % ------------------------------------------------------------------
    e = [ e_f ;      % row 1:     e_i
          e_g ];     % rows 2..:  ebar_i
end


function [e_f, e_g] = fixed_point(F, G, a_all, A_all, ...
                                  w_xi, wbar_xi, w, wbar, R)
% FIXED_POINT  Solve the coupled noise equations (96)-(97) by iteration.
%
%   Inputs
%     F       : 1 x m      noisy function values f_i
%     G       : n x m      noisy gradients g_i
%     a_all   : p x m      column i is a_i
%     A_all   : n x p x m  page i is A_i
%     w_xi    : 1 x m      Taylor-error weights, function values
%     wbar_xi : 1 x m      Taylor-error weights, gradients
%     w       : 1 x m      combined weights  w_i    = w_e_i    + w_xi_i
%     wbar    : 1 x m      combined weights  wbar_i = wbar_e_i + wbar_xi_i
%     R       : p x p      Cholesky factor, S_bar = R' * R
%
%   Outputs
%     e_f : 1 x m      noise terms e_i on the function values
%     e_g : n x m      noise terms ebar_i on the gradients
%
%   Iteration (from the notes):
%     1. Start with e_i = 0 and ebar_i = 0.
%     2. Compute alpha from (S3b), then y = S_bar^{-1} * alpha.
%     3. Update  e_i    = w_xi_i/w_i       * ( f_i - a_i' * y )
%                ebar_i = wbar_xi_i/wbar_i * ( g_i - A_i  * y )
%     4. Repeat from step 2 until the changes are small.

    % ------------------------------------------------------------------
    % Sizes
    % ------------------------------------------------------------------
    n = size(G, 1);            % dimension of the space
    m = size(G, 2);            % number of sample points
    p = size(a_all, 1);        % number of unknowns in y

    % ------------------------------------------------------------------
    % Stopping settings (adjust to taste)
    % ------------------------------------------------------------------
    maxIterations = 1000;      % safety limit
    tolerance     = 1e-12;     % relative change in the noise terms

    % ------------------------------------------------------------------
    % Step 1: start with zero noise
    % ------------------------------------------------------------------
    e_f = zeros(1, m);         % e_i
    e_g = zeros(n, m);         % ebar_i

    converged = false;

    for iteration = 1:maxIterations

        % --------------------------------------------------------------
        % Step 2a: noise-corrected data
        %   r_i    = f_i - e_i
        %   rbar_i = g_i - ebar_i
        % --------------------------------------------------------------
        r    = F - e_f;        % 1 x m
        rbar = G - e_g;        % n x m

        % --------------------------------------------------------------
        % Step 2b: alpha from (S3b)
        %   alpha = sum_i ( w_xi_i * r_i * a_i + wbar_xi_i * A_i' * rbar_i )
        % --------------------------------------------------------------
        alpha = zeros(p, 1);

        for i = 1:m
            a_i = a_all(:, i);         % p x 1
            A_i = A_all(:, :, i);      % n x p

            % Contribution from the function value of point i
            alpha_function = w_xi(i) * r(i) * a_i;

            % Contribution from the gradient of point i
            alpha_gradient = wbar_xi(i) * (A_i' * rbar(:, i));

            alpha = alpha + alpha_function + alpha_gradient;
        end

        % --------------------------------------------------------------
        % Step 2c: y = S_bar^{-1} * alpha, using the factor R
        %   S_bar = R' * R   =>   y = R \ (R' \ alpha)
        % --------------------------------------------------------------
        y = R \ (R' \ alpha);

        % --------------------------------------------------------------
        % Step 3: update the noise terms
        % --------------------------------------------------------------
        e_f_new = zeros(1, m);
        e_g_new = zeros(n, m);

        for i = 1:m
            a_i = a_all(:, i);
            A_i = A_all(:, :, i);

            % --- function noise ---
            % Misfit of the function value, and the share w_xi/w of it
            misfit_f   = F(i) - a_i' * y;           % scalar
            share_f    = w_xi(i) / w(i);            % between 0 and 1
            e_f_new(i) = share_f * misfit_f;

            % --- gradient noise ---
            % Misfit of the gradient, and the share wbar_xi/wbar of it
            misfit_g       = G(:, i) - A_i * y;     % n x 1
            share_g        = wbar_xi(i) / wbar(i);  % between 0 and 1
            e_g_new(:, i)  = share_g * misfit_g;
        end

        % --------------------------------------------------------------
        % Step 4: check how much the noise terms changed
        % --------------------------------------------------------------
        change_f = e_f_new - e_f;
        change_g = e_g_new - e_g;

        % Size of the change and size of the new iterate (all entries)
        changeNorm = norm([change_f; change_g], 'fro');
        newNorm    = norm([e_f_new; e_g_new],   'fro');

        % Accept the new iterate
        e_f = e_f_new;
        e_g = e_g_new;

        % Relative test, with an absolute fallback when the noise is ~0
        if changeNorm <= tolerance * max(newNorm, 1)
            converged = true;
            break;
        end
    end

    % ------------------------------------------------------------------
    % Warn if the iteration hit the limit without converging
    % ------------------------------------------------------------------
    if ~converged
        warning('fixed_point:noConvergence', ...
                ['Fixed point iteration did not converge in %d ', ...
                 'iterations (last change %.3e).'], ...
                maxIterations, changeNorm);
    end
end


function [e_f, e_g] = solve_noise_system(F, G, a_all, A_all, ...
                                         w_xi, wbar_xi, w, wbar, R)
% SOLVE_NOISE_SYSTEM  Solve the stationarity conditions as M*x = b.
%
% The unknown vector is ordered as all function-noise terms first, followed
% by each sample's gradient-noise vector:
%   x = [e_1; ...; e_m; ebar_1; ...; ebar_m].
% The matrix M is the positive-definite Schur complement from the notes.

    n = size(G, 1);
    m = size(G, 2);
    p = size(a_all, 1);
    numberOfUnknowns = m * (n + 1);

    % Build B, where each column contains a sample's contribution to alpha.
    B = zeros(p, numberOfUnknowns);
    beta = zeros(p, 1);
    d = zeros(numberOfUnknowns, 1);
    diagonalW = zeros(numberOfUnknowns, 1);

    for i = 1:m
        a_i = a_all(:, i);
        A_i = A_all(:, :, i);
        gradientColumns = m + (i - 1) * n + (1:n);

        B(:, i) = w_xi(i) * a_i;
        B(:, gradientColumns) = wbar_xi(i) * A_i';

        beta = beta + w_xi(i) * F(i) * a_i ...
            + wbar_xi(i) * (A_i' * G(:, i));

        % d and W use the same function-first, gradient-blocks-second order.
        d(i) = w_xi(i) * F(i);
        d(gradientColumns) = wbar_xi(i) * G(:, i);
        diagonalW(i) = w(i);
        diagonalW(gradientColumns) = wbar(i);
    end

    % Apply S^{-1} using the Cholesky factor S = R' * R; do not form S^{-1}.
    SInvB = R \ (R' \ B);
    SInvBeta = R \ (R' \ beta);

    % Form the Schur-complement system and remove tiny round-off asymmetry.
    M = diag(diagonalW) - B' * SInvB;
    M = 0.5 * (M + M');
    rhs = d - B' * SInvBeta;

    % M is positive definite when the measurement-noise weights are positive.
    % Cholesky gives a stable direct solve without explicitly inverting M.
    Rm = chol(M);
    noise = Rm \ (Rm' \ rhs);

    e_f = noise(1:m)';
    e_g = reshape(noise(m + 1:end), n, m);
end