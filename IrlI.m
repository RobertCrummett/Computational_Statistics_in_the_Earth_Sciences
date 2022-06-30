function [Result] = IrlI(x, Psi, varargin)
% Computes the inverse ilr using the contrast matrix Psi
if nargin == 2
    kappa = 1;
else
    kappa = varargin{1};
end
Result = ClrI(x*Psi, kappa);
end