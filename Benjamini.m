function [RejTh, Index] = Benjamini(pval, varargin)
% Computes the rejection threshold RejTh and the index Index of p-values
% pval that reject the null hypothesis using the Benjamini Hochberg method
% Input Variables:
% pval - vector of p-values (required)
% q - FWER probability level (optional), default is 0.05
% cm - turn on dependance factor (optional), default is one for
% independance

[n p] = size(pval);

switch nargin
    case 1
        q = 0.05;
        cm = 1;
    case 2
        q = varargin{1};
        cm = 1;
    case 3
        q = varargin{1};
        i = 1:max(n, p);
        cm = sum(1./i);
end

if n < p
    pval1 = pval';
else
    pval1 = pval;
end

beta = (1:n)'*q/(cm*n);
pvalo = sort(pval);
i = find(pvalo - beta <= 0);

if isempty(i)
    RejTh = 0;
else
    RejTh = pvalo(i(length(i)));
end

Index = find(pval1 <= RejTh);

end