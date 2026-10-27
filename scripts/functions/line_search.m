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

function print_line_search(lsIter, alpha, fNew, f, dfNew)
    fprintf('LS iteration %2d: alpha = %12d, fvalNew-fval = %12f, dfNew = %12d\n', lsIter, alpha, fNew-f, dfNew);

    if isnan(fNew)
        fprintf('Warning: objective is NaN in line search\n');
    end
end