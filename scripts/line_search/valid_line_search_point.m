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
    finiteOK = all(isfinite(xAll));

    % Kontrollera att samtliga variansparametrar är positiva
    varianceInd = [SigmavInd, SigmawInd, Sigma0Ind];

    varianceValues = xAll(varianceInd);
    badVariance = find(~isfinite(varianceValues) | varianceValues <= 0);

    varianceOK = isempty(badVariance);

    ok = finiteOK && varianceOK;

    % Diagnostik endast när punkten underkänns
    if ~ok
        fprintf('\n[VALIDATION DEBUG] Ogiltig parameterpunkt!\n');

        if ~finiteOK
            badFinite = find(~isfinite(xAll));
            fprintf('Icke-andliga parametrar, globala index:\n');
            disp(badFinite);
        end

        if ~varianceOK
            fprintf('Ogiltiga variansparametrar:\n');

            for j = 1:numel(badVariance)
                k = badVariance(j);
                globalIdx = varianceInd(k);
                value = varianceValues(k);

                if ismember(globalIdx, SigmavInd)
                    group = 'Sigmav';
                    localIdx = find(SigmavInd == globalIdx, 1);
                elseif ismember(globalIdx, SigmawInd)
                    group = 'Sigmaw';
                    localIdx = find(SigmawInd == globalIdx, 1);
                else
                    group = 'Sigma0';
                    localIdx = find(Sigma0Ind == globalIdx, 1);
                end

                fprintf('  %s(%d), global index %d: %.12e\n', ...
                    group, localIdx, globalIdx, value);
            end
        end

        fprintf('[VALIDATION DEBUG] Slut på diagnostik.\n\n');
    end
end
