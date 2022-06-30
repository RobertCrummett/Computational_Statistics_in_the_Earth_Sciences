function [Result] = AInnerProduct(x, y)
% Computes the Aitchison inner product of the compositions x and y
[nx, Dx] = size(x);
[ny, Dy] = size(y);
if Dx ~= Dy
    warning('AInnerProduct: compositions must have the same length')
    return
elseif nx ~= ny
    warning('AInnerProduct: number of observations must be the same')
    return
end

Result = zeros(nx, 1);

x1 = log(x./repmat(geomean(x')', 1, Dx));
y1 = log(y./repmat(geomean(x')', 1, Dy));

for i = 1:nx
    Result(i) = x1(i, :)*y1(i, :)';
end
end