function [pval, lambda] = WilksLambda(h, e, p, n, m)
% Computes the Wilks lambda statistic and double sided p-value from the
% error and hypothesis matrices using Rao (1951) F approximation
%
% Input Variables:
% h - hypothesis matrix
% e - error matrix
% p - dimensionality
% n - hypothesis degrees of freedom
% m - error degrees of freedom

lambda = det(e)/det(e + h);

if p == 1
    f = m*(1 - lambda)/(n*lambda);
    n1 = n;
    n2 = m;
elseif p == 2
    f = (m - 1)*(1 - sqrt(lambda))/(n*sqrt(lambda));
    n1 = 2*n;
    n2 = 2*(m - 1);
elseif n == 1
    f = (n + m - p)*(1 - lambda)/(p*lambda);
    n1 = p;
    n2 = n + m - p;
elseif n == 2
    f = (n + m - p - 1)*(1 - sqrt(lambda))/(p*sqrt(lambda));
    n1 = 2*p;
    n2 = 2*(n + m - p - 1);
else
    n1 = p*n;
    w = n + m - (p + n + 1)/2;
    t = sqrt((p^2*n^2 - 4)/(p^2 + n^2 - 5));
    n2 = w*t - (p*n - 2)/2;
    f = n2*(1 - lambda^(1/t))/(n1*lambda^(1/t));
end

pval = 2*min(fcdf(f, n1, n2), 1 - fcdf(f, n1, n2));

end