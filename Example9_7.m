%% Example 9.7
% Returning to the star data of Example 9.2, use the bounded influence
% estimator to compute a regression fit the the data.

%% Importing Data
clear; close all; clf;
star = importdata('star.dat');

%% Data
y = star(:,2);
x = [star(:,1) ones(size(y))];

%% Chave-Thompson Bounded Influence Algorithm
% Using default leverage point scaling factors.
[tstruct, rstruct, pstruct] = BI(y, x);

%% Truncated Quantile-Quantile Plots
% The bounded influence linear regresion has made some weights zero. This
% means that the parameters have no bearing on the fit, and can thus be
% removed. The truncation of the data means that the resulting q-q plots
% need to be constructed with truncated distributions.

w = rstruct.w(:, 4).*rstruct.what(:, 4);
ii = find(w > 0.00001); %  This picks up data that are nearly eliminated.
res = rstruct.res(:, 4);

[res1, ~] = sort(w(ii).*res(ii));

% The data weighted to zero are assumed evenly distributed between the top
% and the bottom of the distribution.

n = length(w);
m = length(ii);
mr1 = (n - m)/2;
mr2 = n - m - mr1;
i = 1:m;
% Truncated Distriubtion
p = (normcdf((n - mr2 - 1/2)/n) - normcdf((mr1 - 1/2)/n))*(i - 1/2)/m + normcdf((mr1 - 1/2)/n);

% Truncated Residual Q-Q Plot
figure(1);
plot(norminv(p), res1, 'k+', 'LineWidth', 2)

ylabel 'Residual Order Statistic'
xlabel 'Truncated Normal Quantile'
title 'Weighted Residual Q-Q Plot'

% The residual q-q plot is gratifyingly normal.

wx = repmat(w, 1, size(x,2)).*x;
hat = diag(wx*inv(wx'*wx)*wx');
mh1 = 0;
mh2 = n - m;
[hat1, i1] = sort(hat(ii)*sum(w)/2);
% Truncated Distribution
p = betacdf(betainv((n - mh2 - 1/2)/n, size(x,2), size(x,1) - size(x,2)), size(x, 2), size(x, 1) - size(x, 2))*(i - 1/2)/m;

% Truncated Hat Matrix Diagonal Q-Q Plot
figure(2);
plot(betainv(p, size(x,2), size(x,1) - size(x,2)), hat1, 'b+', 'LineWidth', 2)

xlabel 'Truncated Beta (2, 45) Quantile'
ylabel 'Hat Matrix Order Statistic'
title 'Weighted Hat Matrix Diagonal Q-Q Plot'

% The hat matrix diagonal q-q plot is free from leverage points.