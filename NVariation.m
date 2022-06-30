function [T, totvar] = NVariation(x)
% Computes the normalized variation matrix for the composition.
[n, D] = size(x);
k = 1:D;
ratio = zeros(D, D, n);
for i = 1:n
    for j = 1:D
        ratio(j, k, i) = log(x(i, j)./x(i, k));
    end
end
T = var(ratio/sqrt(2), 1, 3);
totvar = sum(sum(T))/D;
end