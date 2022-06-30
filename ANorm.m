function [Result] = ANorm(x)
%Computes the Aitchison norm of the composition x
[n, D] = size(x);
Result = zeros(n, 1);
x1 = log(x./repmat(geomean(x')', 1, D));
for i = 1:n
    Result(i) = sum(x1(i, :).^2);
end
Result = sqrt(Result);
end

