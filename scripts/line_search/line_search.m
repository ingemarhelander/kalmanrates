function [alpha, xNew, fNew, gNew, evals] = line_search(x, f, g, d, K, mu, varInd)
% line_search - Handledarens backtracking-linjesökning (för jämförelse).
%
% Startar med steget från initial_step_length (högst 1, eftersom
% find_search_direction bara garanterar positiva varianser för alpha <= 1)
% och halverar tills F ökar, högst maxLsIter gånger.
%
% evals: alla utvärderade punkter, startpunkten först, en per kolumn
% (evals.X, evals.F, evals.G).

    maxLsIter = 50;
    evals.X = x;  evals.F = f;  evals.G = g;
    alpha = initial_step_length(g, d);
    xNew = x + alpha*d;
    [fNew, gNew] = kalmanObj(xNew, K, mu, varInd);
    evals.X(:, end+1) = xNew;  evals.F(1, end+1) = fNew;  evals.G(:, end+1) = gNew;
    lsIter = 1;
    while true
        print_line_search(lsIter, alpha, fNew, f, gNew'*d)
        if fNew > f
            break;                   % förbättring: klart
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
