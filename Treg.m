function [b, bse, t, pval] = Treg(y, x, w)
% Computes the t table for a linear regression.
% 
% The null hypothesis is rejected if the p-value for the test statistic is
% below the Type 1 error level alpha.
%
% It is important to note that the test if for betahat_j = 0 in a linear
% model that includes all the remianing predictors. So it is not an absolute
% test of betahat_j. 

p = size(x,2);
i = find(w ~= 0);
b = (repmat(w(i),1,p).*x(i,:))\(w(i).*y(i));

yhat = x(i,:)*b;
n = length(i);

% Sample Statistics
zeta = inv((repmat(w(i), 1, p)'.*x(i,:)')*(repmat(w(i), 1, p).*x(i,:)));
sigma2 = (w(i)'.*(y(i) - yhat)')*(w(i).*(y(i) - yhat))/(n-p);

bse = sqrt(sigma2*diag(zeta));
t = b./bse;

% P-Value
pval = 2*(1 - tcdf(abs(t), n - p));

end