function ok = valid_line_search_point(xNew, K, varInd)

    % Återskapa hela parametervektorn på samma sätt som i kalmanObj
    thetaFInd = 1:length(K.thetaF);
    SigmavInd = (1:length(K.Sigmav)) + thetaFInd(end);
    SigmawInd = (1:length(K.Sigmaw)) + SigmavInd(end);
    axInd = (1:length(K.ax)) + SigmawInd(end);
    Sigma0Ind = (1:length(K.Sigma0)) + axInd(end);

    xAll = [K.thetaF; K.Sigmav; K.Sigmaw; ...
            K.ax; K.Sigma0; K.Ef(:)];

    xAll(varInd) = xNew;

    % Kontrollera att alla parametrar är ändliga
    ok = all(isfinite(xAll));

    % Kontrollera att samtliga variansparametrar är positiva
    varianceInd = [SigmavInd, SigmawInd, Sigma0Ind];

    ok = ok && all(xAll(varianceInd) > 0);
end