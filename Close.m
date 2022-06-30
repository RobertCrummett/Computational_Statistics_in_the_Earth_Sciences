function [Result] = Close(x, varargin)
% Closes a composition x to the constant kappa that defaults to one if not
% provided
if nargin == 1
    kappa = 1;
else
    kappa = varargin{1};
end
[~, D] = size(x);
Result = kappa*x./(repmat(sum(x')', 1, D));
end
