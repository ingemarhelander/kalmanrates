function [f, grad] = kalmanObjMin(x, K, mu, varInd)

[f, grad] = kalmanObj(x, K, mu, varInd);
f = -f;
grad = -grad;

