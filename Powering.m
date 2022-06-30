function [Result] = Powering(alpha, x, varargin)
% Computes the power transformation of a compositional vector x
if nargin == 2
    kappa = 1;
else
    kappa = varargin{1};
end
Result = Close(x.^alpha, kappa);
end