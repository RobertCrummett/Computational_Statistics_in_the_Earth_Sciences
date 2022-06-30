function [Result] = ADistance(x, y)
% Computes the Aitchison distance between compositions x and y
Result = ANorm(Perturbation(x, 1./y));
end