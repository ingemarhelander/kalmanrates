function [alpha, xNew, fNew, gNew, evals] = line_search(x, f, g, d, K, mu, varInd)
% Backtracking line search: start with a (capped) step, halve until the
% objective improves or the maximum number of halvings is reached.
% evals holds every point evaluated in this step (start point included),
% one per column: evals.X, evals.F, evals.G.
    maxLsIter = 50;
    evals.X = x;  evals.F = f;  evals.G = g;
    alpha = first_step_length(g, d);
    xNew = x + alpha*d;
    [fNew, gNew] = kalmanObj(xNew, K, mu, varInd);
    evals.X(:, end+1) = xNew;  evals.F(1, end+1) = fNew;  evals.G(:, end+1) = gNew;
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
        evals.X(:, end+1) = xNew;  evals.F(1, end+1) = fNew;  evals.G(:, end+1) = gNew;
        lsIter = lsIter + 1;
        if lsIter > maxLsIter
            break;
        end
    end
end

function print_line_search(lsIter, alpha, fNew, f, dfNew)
    fprintf('LS iteration %2d: alpha = %12d, fvalNew-fval = %12f, dfNew = %12d\n', lsIter, alpha, fNew-f, dfNew);

    if isnan(fNew)
        fprintf('Warning: objective is NaN in line search\n');
    end
end