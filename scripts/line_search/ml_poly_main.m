function [q, fit, info] = ml_poly_main(alpha, F, S, qPrev, method)
% ml_poly_main - ML-polynom längs linjen (ML-boxen i line search).
%
% Ingång till ML-delen. Väljer polynomgrad med en av fyra metoder:
%   'score'     mult_q_poly_fits + select_poly_degree, viktad poäng
%               (vikterna sätts i select_poly_degree)
%   'threshold' mult_q_poly_fits + select_poly_degree, lägsta grad som
%               klarar alla tröskelvärden
%   'lr'        select_poly_degree_lr, likelihoodkvottest
%   'lag1'      select_poly_degree_lag1, lag-1-autokorrelation
%
%   alpha   steglängder (bara punkter där F och S är ändliga)
%   F       funktionsvärden
%   S       lutning längs linjen, g'*d
%   qPrev   grad från förra varvet, eller [] första gången (valfri, bara 'lag1')
%   method  se ovan (valfri, standard 'lr')
%   q       vald grad, fit är ML-anpassningen för den graden
%   info    (valfri) detaljer per provad grad, se respektive funktion

    if nargin < 4
        qPrev = [];
    end
    if nargin < 5 || isempty(method)
        method = 'lr';
    end

    % info begärs bara när den behövs: select_poly_degree_lag1 räknar då tau
    % även för sista graden, vilket annars är onödigt arbete.
    info = [];
    switch method
        case {'score', 'threshold'}
            [q, fit, info] = select_with_diagnostics(alpha, F, S, method);
        case 'lag1'
            if nargout > 2
                [q, fit, info] = select_poly_degree_lag1(alpha, F, S, qPrev);
            else
                [q, fit] = select_poly_degree_lag1(alpha, F, S, qPrev);
            end
        case 'lr'
            [q, fit, info] = select_poly_degree_lr(alpha, F, S);
        otherwise
            error(['ml_poly_main: okänd metod ''%s'', använd ''score'', ''threshold'', ' ...
                   '''lr'' eller ''lag1''.'], method);
    end
end

function [q, fit, info] = select_with_diagnostics(alpha, F, S, method)
% Anpassar grad 1..qMax med mult_q_poly_fits och låter select_poly_degree
% välja. Högsta grad är som i select_poly_degree_lr, m - 2 (minst 2, högst
% 6), så att varken F eller S kan anpassas exakt.
%
% Reserv q = 2 när
%   - färre än 3 punkter finns: trendtestet i diagnostiken kräver minst 3
%     punkter, så ingen grad blir giltig,
%   - select_poly_degree inte väljer någon grad ('threshold' när ingen grad
%     klarar trösklarna). info.status anger vilket.
    m = numel(alpha);
    qMax = max(2, min(6, m - 2));
    fits = mult_q_poly_fits(alpha, F, S, qMax, false);

    if m < 3
        q = 2;
        fit = fits{2};
        info = struct('method', method, 'status', 'too_few_points_fallback_q2');
        return
    end
    % Grad 1 har inget maximum och kan inte användas av find_step_length
    [q, fit, info] = select_poly_degree(fits, method, struct('minDegree', 2));
    if isempty(q)
        q = 2;
        fit = fits{2};
        info.status = [info.status '_fallback_q2'];
    end
end
