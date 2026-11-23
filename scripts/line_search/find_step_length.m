function alphaNew = find_step_length(fit, h, alphaMax)
% find_step_length - Nästa steglängd: polynomets maximum.
%
% Evaluerar det anpassade polynomet i 50 jämnt fördelade punkter
% inom intervallet [0, alphaMax].
%
% Punkten med högst predikterat funktionsvärde returneras.
% Metoden hittar därmed även maximum när polynomet initialt har
% negativ lutning.
%
% Ingen kalmanObj körs här.
%
% In:   fit       ML-anpassningen (från select_poly_degree)
%       h         steglängd (behålls för kompatibilitet)
%       alphaMax  längsta tillåtna steg
% Ut:   alphaNew  steglängden med högst predikterat funktionsvärde

    % Skapa 50 jämnt fördelade testpunkter
    steps = linspace(0, alphaMax, 50);

    % Evaluera polynomet i samtliga punkter
    values = poly_value(fit, steps);

    % Välj punkten med högst predikterat funktionsvärde
    [~, k] = max(values);

    alphaNew = steps(k);

end
