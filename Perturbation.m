function [Result] = Perturbation(x, y, varargin)
% Performs  the compositional operation perturbation on two compositions
if nargin == 2
    kappa = 1;
else
    kappa = varargin{1};
end

[nx, Dx] = size(x);
[ny, Dy] = size(y);
x1 = x;
y1 = y;

if Dy == 1
    y1 = repmat(y, 1, Dx);
elseif nx == 1
    x1 = repmat(x, ny, 1);
elseif ny == 1
    y1 = repmat(y, nx, 1);
end

if Dx ~= Dy
    warning('Perturbation: compositions must have the same length');
    return
elseif nx ~= ny
    warning('Perturbation: number of observations must be the same');
    return
end

Result = Close(x1.*y1, kappa);
end