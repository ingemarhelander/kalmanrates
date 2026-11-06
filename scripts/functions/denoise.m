function d = denoise(x, f, g, D, evalHist, noisy_size, K, mu, varInd, xLbd)
% Placeholder: returns x, f, g unchanged.
    columnNorms = sqrt(sum(abs(D).^2, 1));

    %Check that no columns are of length 0.
    if any(columnNorms == 0)
        error('denoise:ZeroColumnInD', 'D cannot contain zero columns.');
    end

    % Normalize the columns of D to get an orthonormal basis for the span of the past search directions.
    D = D ./ columnNorms;

    fprintf('\nNoisy optimization: last %d search directions \n', noisy_size);
    %fprintf('denoise: evalhist has %d steps, evaluations per step = [%s], total = %d \n', ...
    %    numel(evalHist), num2str(cellfun(@(e) size(e.X,2), evalHist)), sum(cellfun(@(e) size(e.X,2), evalHist)))

    % Get the X, F andf G values from the evalHist.
    X = [];
    F = [];
    G = [];
    for k = 1:numel(evalHist)
        X = [X, evalHist{k}.X];
        F = [F, evalHist{k}.F];
        G = [G, evalHist{k}.G];
    end
    % Eval j corresponds to X(:,j), F(j) and G(:,j)

    fprintf('[denoise] size(D)=%s, size(X)=%s, size(F)=%s, size(G)=%s\n', ...
        mat2str(size(D)), mat2str(size(X)), mat2str(size(F)), mat2str(size(G)));

    [Q, q, c, e] = get_quadratic_approx_in_subspace(X, F, G, D);

    % Check that q is finite.
    if any(~isfinite(q))
        error('denoise:NonFiniteDirection', 'The quadratic model produced a non-finite direction.');
    end

    % Compute the search direction d from the quadratic model.
    qNorm = norm(q);
    if qNorm == 0
        d = zeros(size(q));
        return;
    end
    d = q / qNorm;

    % Make sure tyhe search direction d is feasible with respect to the lower bounds xLbd.
    alphaLbd = ones(size(x)) * inf;
    hitsBound = (d < 0) & (xLbd ~= -inf);
    alphaLbd(hitsBound) = -(x(hitsBound) - xLbd(hitsBound)) ./ d(hitsBound);
    alphaLbdMin = min(alphaLbd);
    if alphaLbdMin < 1
        d = d * alphaLbdMin * 0.99;
    end
end


function [Q, q, c, e] = get_quadratic_approx_in_subspace(X, F, G, D)
% Project X and G onto the span of D to get a smaller problem.
    Z = D' * X;
    Gsmall = D' * G;
    fprintf('[denoise] approximation inputs: size(Z)=%s, size(F)=%s, size(Gsmall)=%s\n', ...
        mat2str(size(Z)), mat2str(size(F)), mat2str(size(Gsmall)));

    %Call the function that computes the quadratic approximation.
    [Qsmall, qsmall, c, esmall] = compute_quadratic_approximation(Z, F, Gsmall);


    fprintf('[denoise] approximation outputs: size(Qsmall)=%s, size(qsmall)=%s, size(c)=%s, size(esmall)=%s\n', ...
        mat2str(size(Qsmall)), mat2str(size(qsmall)), mat2str(size(c)), mat2str(size(esmall)));

    % Transform the small quadratic model back to the original space.
    Q = D * Qsmall * D';
    q = D * qsmall;
    e = [esmall(1, :); D * esmall(2:end, :)];
end


function [Q, q, c, e] = compute_quadratic_approximation(X, F, G)
% COMPUTE_QUADRATIC_APPROXIMATION  Fit one quadratic model to noisy f and g.
%
%   Inputs
%     X : n x m   the m sample points (one per column)
%     F : 1 x m   noisy function values f_i
%     G : n x m   noisy gradients g_i (one per column)
%
%   Outputs
%     Q : n x n   Hessian of the model
%     q : n x 1   gradient of the model at xbar
%     c : scalar  value of the model at xbar
%     e : (1+n) x m  noise terms. Row 1 is e_i (function noise),
%                    rows 2..n+1 are ebar_i (gradient noise).
%
%   The model is centred at xbar = last column of X, and
%       f_i = 1/2 dx_i' Q dx_i + q' dx_i + c + e_i + xi_i
%       g_i = Q dx_i + q + ebar_i + xibar_i

    % ------------------------------------------------------------------
    % 0. Sizes and shifted points
    % ------------------------------------------------------------------
    n = size(X, 1);            % dimension of the space
    m = size(X, 2);            % number of sample points

    % Shift all points so that xbar (the last point) is the origin.
    deltaX = X - X(:, end);    % column i is dx_i = x_i - xbar

    % Number of unique elements in the symmetric matrix Q.
    nVech = n * (n + 1) / 2;

    % Total number of unknowns in y = [vech(Q); q; c].
    p = nVech + n + 1;

    % ------------------------------------------------------------------
    % 0b. Weights (PLACEHOLDERS: the notes do not say how to choose them)
    % ------------------------------------------------------------------
    % w_xi(i)     : penalty on the Taylor error xi_i in the function value
    % wbar_xi(i)  : penalty on the Taylor error xibar_i in the gradient
    % w_e(i)      : penalty on the noise e_i in the function value
    % wbar_e(i)   : penalty on the noise ebar_i in the gradient
    w_xi    = ones(1, m);
    wbar_xi = ones(1, m);
    w_e     = ones(1, m);
    wbar_e  = ones(1, m);

    % Combined weights from Step 7:  w_i = w_e_i + w_xi_i
    w    = w_e    + w_xi;
    wbar = wbar_e + wbar_xi;

    % ------------------------------------------------------------------
    % Step 1: Compact notation
    % ------------------------------------------------------------------

    % --- 1a. Duplication matrix D, defined by vec(Q) = D * vech(Q) -----
    % vech stacks the lower triangle of Q column by column.
    D = zeros(n * n, nVech);

    k = 0;                                % running index into vech(Q)
    for col = 1:n
        for row = col:n                   % lower triangle: row >= col
            k = k + 1;

            % Position of Q(row, col) inside vec(Q) (column-major)
            idxLower = (col - 1) * n + row;

            % Position of the mirrored element Q(col, row) inside vec(Q)
            idxUpper = (row - 1) * n + col;

            % Both positions in vec(Q) take their value from vech(Q)(k)
            D(idxLower, k) = 1;
            D(idxUpper, k) = 1;
        end
    end

    % --- 1b. Build a_i (p x 1) and A_i (n x p) for every point ---------
    % a_i  such that  1/2 dx' Q dx + q' dx + c  =  a_i' * y
    % A_i  such that  Q dx + q                  =  A_i * y
    a_all = zeros(p, m);         % column i is a_i
    A_all = zeros(n, p, m);      % page i is A_i

    for i = 1:m

        dx = deltaX(:, i);       % dx_i as a column vector

        % ---- a_i = [ D' * vec(1/2 dx dx') ; dx ; 1 ] ----
        outerProduct = dx * dx';                  % n x n matrix dx dx'
        halfOuter    = 0.5 * outerProduct;        % 1/2 dx dx'
        halfOuterVec = halfOuter(:);              % vec(.), n^2 x 1
        topPart      = D' * halfOuterVec;         % nVech x 1

        a_i = [ topPart ;
                dx      ;
                1       ];

        a_all(:, i) = a_i;

        % ---- A_i = [ Atilde_i * D , I , 0 ] ----
        % Q*dx = kron(dx', I_n) * vec(Q), so Atilde_i = kron(dx', I_n).
        Atilde_i = kron(dx', eye(n));             % n x n^2

        A_i = [ Atilde_i * D , ...                % n x nVech
                eye(n)       , ...                % n x n  (picks out q)
                zeros(n, 1) ];                    % n x 1  (c does not enter)

        A_all(:, :, i) = A_i;
    end

    % ------------------------------------------------------------------
    % Step 3 (S3c): S_bar does not depend on the noise, so build it once
    % ------------------------------------------------------------------
    % S_bar = sum_i ( w_xi_i * a_i a_i' + wbar_xi_i * A_i' A_i )
    S_bar = zeros(p, p);

    for i = 1:m
        a_i = a_all(:, i);
        A_i = A_all(:, :, i);

        % Contribution from the function value of point i
        S_function = w_xi(i) * (a_i * a_i');

        % Contribution from the gradient of point i
        S_gradient = wbar_xi(i) * (A_i' * A_i);

        S_bar = S_bar + S_function + S_gradient;
    end

    % Make sure S_bar is exactly symmetric (removes rounding asymmetry)
    S_bar = 0.5 * (S_bar + S_bar');

    % Factorise once. S_bar is constant, so the factor can be reused in
    % every iteration. Needs S_bar positive definite, which requires
    % enough points in general position (see note below).
    %   S_bar = R' * R,  so  S_bar \ b = R \ (R' \ b)
    R = chol(S_bar);

    % ------------------------------------------------------------------
    % Steps 5-7: solve the coupled equations (96)-(97) for the noise
    % ------------------------------------------------------------------
    % Returns e_f (1 x m) and e_g (n x m).
    %[e_f, e_g] = fixed_point(F, G, ...
    %                         a_all, A_all, ...
    %                         w_xi, wbar_xi, w, wbar, ...
    %                         R);
    [e_f, e_g] = solve_noise_system(F, G, ...
                                a_all, A_all, ...
                                w_xi, wbar_xi, w, wbar, ...
                                R);

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

    % ------------------------------------------------------------------
    % Unpack y = [vech(Q); q; c] into Q, q, c
    % ------------------------------------------------------------------
    vechQ = y(1 : nVech);
    q     = y(nVech + 1 : nVech + n);
    c     = y(end);

    % vec(Q) = D * vech(Q), then reshape the vector into an n x n matrix
    vecQ = D * vechQ;
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
    rhs = -(d + B' * SInvBeta);

    % M is positive definite when the measurement-noise weights are positive.
    % Cholesky gives a stable direct solve without explicitly inverting M.
    Rm = chol(M);
    noise = Rm \ (Rm' \ rhs);

    e_f = noise(1:m)';
    e_g = reshape(noise(m + 1:end), n, m);
end