function config = modelRunConfig()
% Shared data limit for initialization and saved-model optimization.
% Use Inf on the full-data computer, or a positive integer for laptop tests.
config.maxObservations = 250;
end
