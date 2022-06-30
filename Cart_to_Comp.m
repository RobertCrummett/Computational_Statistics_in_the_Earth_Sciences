function [Result] = Cart_to_Comp(xx)
%Converts Cartesian components xx on the ternary with unit edges to a
%composition Result
n = size(xx, 1);
Result = zeros(n, 3);
for i = 1:n
    Result(i, 1:3) = [2*xx(i, 2)/sqrt(3) xx(i, 1) - xx(i, 2)/sqrt(3) ...
                      1 - xx(i, 1) - xx(i, 2)/sqrt(3)];
end

