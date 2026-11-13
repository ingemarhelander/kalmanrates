function [q, fit, info] = ml_poly_main(alpha, F, S, qPrev, method)
% ml_poly_main - ML-polynom längs linjen (ML-boxen i line search).
%
% Ingång till ML-delen. Väljer polynomgrad med en av två metoder:
%   'lr'    select_poly_degree_lr, likelihoodkvottest (standard)
%   'lag1'  select_poly_degree_lag1, lag-1-autokorrelation
%
%   alpha   steglängder (bara punkter där F och S är ändliga)
%   F       funktionsvärden
%   S       lutning längs linjen, g'*d
%   qPrev   grad från förra varvet, eller [] första gången (valfri, bara 'lag1')
%   method  'lr' eller 'lag1' (valfri, standard 'lr')
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
        case 'lag1'
            if nargout > 2
                [q, fit, info] = select_poly_degree_lag1(alpha, F, S, qPrev);
            else
                [q, fit] = select_poly_degree_lag1(alpha, F, S, qPrev);
            end
        case 'lr'
            [q, fit, info] = select_poly_degree_lr(alpha, F, S);
        otherwise
            error('ml_poly_main: okänd metod ''%s'', använd ''lr'' eller ''lag1''.', method);
    end
end
