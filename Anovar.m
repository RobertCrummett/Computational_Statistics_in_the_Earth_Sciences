function [a] = Anovar(y, x, w)
% Computes the ANOVA table for a regression problem.
%
% Variables:
% y - response variable
% x - predictor variable
% w -- vector of weights lying between 0 and 1
% one - set to zero or one as x does not/does contain a collumn of ones
%
% Null hypothesis: Beta = 0
% When the null hypothesis is true, the elements of b should be distributed
% according to F(p-1,n-p).
% The test is quite blunt because one is usually interested in testing a
% certain paramter. This effectively tests the applicablility of regression
% in general. 
%
% ANOVA table (a):
% Rows: [Due to regression (explained)],[Due to error (residual)],[Total]
% Columns: [dof],[Sum of Squares],[Mean Squares],[F],[p-Value]

p = size(x, 2);
one = 0;

for i = p:-1:1
    ii = find(x(:,i) ~= 1);
    if isempty(ii)
        one = 1;
        break
    end
end

i = find(w ~= 0);
b = (repmat(w(i), 1, p).*x(i,:))\(w(i).*y(i));
yhat = x(i,:)*b;

n = length(i);
wy = w(i).*y(i);
wyhat = w(i).*yhat;

% Total Sum of Squares
tss = (wy - mean(wy))'*(wy - mean(wy));
% Residual (Between) Sum of Squares
rss = (wy - wyhat)'*(wy - wyhat);
% Explained (Within) Sum of Squares
ess = wyhat'*wyhat - n*mean(wy)^2;

me = ess/(p - one);
mr = rss/(n - p);
f = me/mr;
% P-Value
fpval = 1 - fcdf(f, p - one, n - p);

% ANOVA table
a = [p - one ess me f fpval;
    n - p rss mr 0 0;
    n - one tss 0 0 0];

end