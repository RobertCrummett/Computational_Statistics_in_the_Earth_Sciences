function [Result] = Comp_to_Cart(comp)
%Converts a closed composition comp to Cartesian components Result on a
%ternary with unit edges
[n D] = size(comp);
if D ~= 3
    warning('Comp_to_Cart: composition must have 3 parts')
    return
end
Result = zeros(n, 2);
for i = 1:n
    x = (comp(i, 1) + 2*comp(i, 2))/(2*(comp(i, 1) + comp(i, 2) + ...
        comp(i, 3)));
    y = sqrt(3)*comp(i, 1)/(2*(comp(i, 1) + comp(i, 2) + comp(i, 3)));
    Result(i,:) = [x y];
end


