function fit = poly_fit_ml(alpha, F, S, q)
% poly_fit_ml - ML-anpassning av ett polynom längs linjen (Req. 23).
%
% TILLFÄLLIG VERSION: vanlig minsta kvadrat utan vikter, bara för att
% gradvalet ska gå att köra och testa. Oskar ersätter den med ML-skattningen.
%
% In (kolumner med m element, sorterade efter alpha):
%   alpha  steglängder
%   F      funktionsvärden
%   S      lutning längs linjen, g'*d (derivata med avseende på alpha)
%   q      polynomgrad
%
% Ut:
%   fit.y           polynomets koefficienter (format valfritt, används av maxsökningen)
%   fit.rF, fit.rG  residualer för F och S (m x 1), används av select_poly_degree
%   fit.wF, fit.wG  vikter 1/sigmaF^2 och 1/sigmaG^2, används av select_poly_degree

    a = alpha.^(0:q);                                   % rader [1, α, α², ...]
    b = [zeros(size(alpha)), (1:q).*alpha.^(0:q-1)];    % derivatan av raderna i a
    fit.y  = [a; b] \ [F; S];
    fit.rF = F - a*fit.y;
    fit.rG = S - b*fit.y;
    fit.wF = 1;
    fit.wG = 1;
end
