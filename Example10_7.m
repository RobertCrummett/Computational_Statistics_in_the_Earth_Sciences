%% Example 10.7
% The last two variables will be regressed on the first nine, and nackward
% stepwise selection will be applied to cull the set of predictors.

%% Importing Data
clear; clf; close all;
the = importdata('the.dat');

%% Data
x = [the(:, 1:9) ones(size(the(:,1)))];
y = the(:, 10:11);
[n, m] = size(x);

%% Multivariate Regression
beta = x\y;

e = (y - x*beta)'*y;
h = beta'*x'*y - length(y)*mean(y)'*mean(y);

[pval, lambda] = WilksLambda(h, e, 2, m - 1, n - m - 2);

% The p-value is very small. Hence the null hypothesis that all the
% regression coefficients are zero is rejected.

%% Backward Elimination
% The backward eliminaton procedure will be applied to determine if all the
% predictor variables are required.

cv = 0.05; % p-value to achieve below
stop = 0;
removed = [];

while ~stop
    q = 2;
    nh = 1;
    p = size(x,2) - 1;
    ne = length(y) - p - 1;
    lambda2 = [];

    for i = 1:p
        if i == 1
            x1 = x(:, 2:p + 1);
        elseif i < p
            x1 = [x(:, 1:i) x(:, (i + 2):(p + 1))];
        else
            x1 = [x(:, 1:(p - 1)) x(:, p + 1)];
        end

        beta1 = x1\y;

        e = (y - x1*beta1)'*y;
        h = beta1'*x1'*y - length(y)*mean(y)'*mean(y);

        lambda1 = lambda*det(e + h)/det(e);
        lambda2 = [lambda2 lambda1];
    end

    f = (ne - p + 1)*(1 - lambda2)./(p*lambda2);

    pval = 2*min(fcdf(f, p, ne - p + 1), 1 - fcdf(f, p, ne - p + 1));

    [~, idx] = max(lambda2);
    x(:, idx) = [];

    if isempty(removed)
        removed = [removed idx];
    elseif idx >= removed(end)
        removed = [removed (idx + 1)];
    elseif idx < removed(end)
        removed = [removed idx];
    end

    if pval < cv*ones(size(pval))
        break
    end

end

% Now the p-values are all below 0.05, and the process terminates.
% Four different predictors were eliminated. 