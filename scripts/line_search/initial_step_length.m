function alpha = first_step_length(g, d)
% first_step_length - Step length to the first new point of the line search.
%
% The line search starts with one evaluated point (alpha = 0). It needs one
% more point along d before a polynomial can be fitted (Req. 22).

% Choose alpha so that the tangent L(alpha) = f + alpha*(g'*d) predicts
% the increase deltaF. A steep slope gives a short step, a flat slope a long step.

% The result is limited to [alphaMin, 1]. Steps above 1 are not allowed since
% find_search_direction only guarantees the lower bounds for alpha <= 1.

    deltaF   = 1E2;   % Predicted increase in the objective value at the first new point of the line search
    alphaMin = 1E-5;  % Shortest allowed first step 

    slope = g'*d;     % Directional derivative at the start point

    % If slope < 0 the division gives a negative step, which alphaMin replaces.
    % If slope = 0 the division gives Inf, which the upper limit 1 replaces.
    
    alpha = deltaF/slope;
    
    alpha = min(alpha, 1);
    alpha = max(alpha, alphaMin);
end
