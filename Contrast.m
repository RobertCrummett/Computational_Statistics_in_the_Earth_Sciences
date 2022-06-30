function [Result] = Contrast(sbp)
% Computes the contrast matrix from the sequential binary partition sbp
[Dm1, D] = size(sbp);
if D ~= Dm1 + 1
    warning('Contrast: sbp has wrong number of rows')
    return
end
Result = zeros(Dm1, D);
for i = 1:Dm1
    r = find(sbp(i,:) == 1);
    rs = sum(sbp(i, r));
    s = find(sbp(i, :) == -1);
    ss = -sum(sbp(i, s));
    psip = sqrt(ss/(rs*(rs + ss)));
    psim = -sqrt(rs/(ss*(rs + ss)));
    Result(i,r) = psip;
    Result(i,s) = psim;
end
end