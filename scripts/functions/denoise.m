function [x, f, g] = denoise(x, f, g, D, evalHist, noisy_size, K, mu, varInd, xLbd)
% Placeholder: returns x, f, g unchanged.
    fprintf('\nNoisy optimization: last %d search directions \n', noisy_size);
    %fprintf('denoise: evalhist has %d steps, evaluations per step = [%s], total = %d \n', ...
    %    numel(evalHist), num2str(cellfun(@(e) size(e.X,2), evalHist)), sum(cellfun(@(e) size(e.X,2), evalHist)))
    X = [];
    F = [];
    G = [];
    for k = 1:numel(evalHist)
        X = [X, evalHist{k}.X];
        F = [F, evalHist{k}.F];
        G = [G, evalHist{k}.G];
    end
    %Eval j corresponds to X(:,j), F(j) and G(:,j)
end