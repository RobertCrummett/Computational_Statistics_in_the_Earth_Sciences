function [Result] = ClrI(x, varargin)
% Computes the inverse clr of x
if nargin == 1
    kappa = 1;
else
    kappa = varargin{1};
end
Result = Close(exp(x), kappa);
end