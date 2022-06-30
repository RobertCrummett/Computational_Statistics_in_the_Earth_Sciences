function [Result] = Clr(x)
% Applies the centered log ratio transformation to the composition x
[~, D] = size(x);
Result = log(x./repmat(geomean(x')', 1, D));
end