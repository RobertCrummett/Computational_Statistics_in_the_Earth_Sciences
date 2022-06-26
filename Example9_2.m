%% Example 9.2
% A least squares line is fit to the data after adding a column of ones to
% the predictor matrix to allow for a nonzero intercept, and then the
% result is plotted.

%% Importing Data
clear; close all; clf;
star = importdata('star.dat');

%% Data
y = star(:,2);
x = [star(:,1) ones(size(y))];

%% Linear Regression
b = x\y; % coefficent estimate
x1 = 3.0 : 0.1 : 5.0;
y1 = b(1)*x1 + b(2);

%% Plotting Least Squares Fit
figure(1);
plot(x(:,1), y, 'k+', 'LineWidth',1.5)
hold on
plot(x1, y1, 'r--','LineWidth',2)
hold off

xlabel 'Log Temperature';
ylabel 'Log Light Intensity'
title 'Star Data and Best-Fit Line Using All the Data'

%% Quality Check
% Runs Test
r = y - x*b; % residuals
[h, p] = runstest(r);

% D-W Test
pval = dwtest(r, x);

% The Durbin-Watson test shows that the data are randomand uncorrelated,
% which is counter intuitive because the line that should best fit the data
% should have a positive slope.

%% ANOVA and t-Statistic Analysis
w = ones(size(y));

% ANOVA
a = Anovar(y, x, w);

% t-Statistic
[b, bse, t, tprob] = Treg(y, x, w);

% Correlation Coefficient
r2 = a(1,2)/a(3,2);

% The ANOVA p-value rejects the need for regression. The regression fits
% nothing very well.

%% Q-Q Plots Examined
i = 1:length(y);
p = (i - 1/2)/length(y);
[r1, i1] = sort(r);

figure(2);
plot(norminv(p), r1, '+k', 'LineWidth', 2)

xlabel 'Normal Qunatiles'
ylabel 'Residual Order Statistic'
title 'Q-Q Plot for the Star Residuals Against Gaussian Quantiles'

% The q-q plot is approximately straight and there are no outliers,
% although there is weak systematic curvature.

h = x*inv(x'*x)*x';
hii = diag(h);
[h1, i1] = sort(hii);

figure(3);
plot(betainv(p, 2, length(y) - 2), h1, 'b+', 'LineWidth', 2)

xlabel 'Beta(2, 45) Quantiles'
ylabel 'Hat Matrix Order Statistic'
title(["Q-Q Plot for the Hat Matrix Diagonal Normalized"; ...
    "by Expected Value Against Beta(2, 45) Quantiles"])

% The hat matrix shows the presence of four high -leverage points and
% possibly a fifth. These can be removed and the regression proceedure can
% be redone.

%% Removing Data
points2remove = 4;
w(i1(end - points2remove - 1 : end), :) = 0; % accomplished by weighting them as 0

%% Regression on Truncated Data
% ANOVA
a = Anovar(y, x, w);

% t-Statistic Test
[b, bse, t, tprob] = Treg(y, x, w);
r2 = a(1,2)/a(3,2); % correlation coefficent

% The test concludes that both coefficients are significant.

% Runs Test
i = find(w ~= 0);
r = y(i) - x(i,:)*b; % residuals
[h, p] = runstest(r);

% D-W Test
x1 = x(i, :);
p = dwtest(r, x1);

% The Durbin-Watson test rejects the null hypothesis that the data are 
% uncorrelated, while the runs test accepts the null hypothesis that the 
% residuals are random.

%% Linear Regression for Truncated Data
x1 = 3.0 : 0.1 : 5.0;
y1 = b(1)*x1 + b(2);

%% Plotting Least Squares Fit for Truncated Data
figure(4);
plot(x(i,1), y(i), 'k+', 'LineWidth',1.5)
hold on
plot(x1, y1, 'r--','LineWidth',2)
hold off

xlabel 'Log Temperature';
ylabel 'Log Light Intensity'
title 'Star Data and Best-Fit Line Using Truncated Data'
xlim([3.5 5]); 
ylim([3 6.5]);

%% Q-Q Plots for the Truncated Data
% The q-q plots of the regression residuals and hat matrix may be examined
% using the appropriate truncateed forms of the distributions.

n = length(y);
m = length(i);
ii = 1 : m;

% Assume 2 values truncated from bottom and two values from the top of the
% normal distirbution.

mr1 = 2;
mr2 = 2;
p = (normcdf(norminv((n - mr2 - 1/2)/n)) - normcdf(norminv((mr1 - 1/2)/n)))*...
    (ii - 1/2)/m + normcdf(norminv((mr1 - 1/2)/n));
[r1, i1] = sort(r);

figure(5);
plot(norminv(p), r1, 'k+', 'LineWidth', 2)
hold on
[~,~,~,cv] = kstest(r1);
plot(norminv(p), norminv(p + cv, mean(r1), std(r1)), 'r--', 'LineWidth', 2)
plot(norminv(p), norminv(p - cv, mean(r1), std(r1)), 'r--', 'LineWidth', 2)
hold off

xlabel 'Truncated Normal Qunatiles'
ylabel 'Residual Order Statistic'
legend 'Data' 'C_{\alpha} = 0.95' 'Location' 'se'
title 'Q-Q Plot for the Star Residuals Against Truncated Gaussian Quantiles'

h = x(i,:)*inv(x(i,:)'*x(i,:))*x(i,:)';
hii = m/2*diag(h);
[h1, i1] = sort(hii);

% Truncation of the beta distribution is only at the top.

mh2 = 4;
p = betacdf(betainv((n - mh2 - 1/2)/n, 2, n - 2), 2, n - 2)*(ii - 1/2)/m;

figure(6);
plot(betainv(p, 2, n - 2), h1, 'b+', 'LineWidth', 2)

xlabel 'Truncated Beta(2, 45) Quantiles'
ylabel 'Hat Matrix Order Statistic'
title(["Q-Q Plot for the Hat Matrix Diagonal Normalized"; ...
    "by Expected Value Against Truncated Beta(2, 45) Quantiles"])

%% Scatter Plot of Final Residuals
figure(7);
scatter(i, r, '+', 'LineWidth', 1.5)

ylim([-1 1.5])
xlabel 'Data Index'
ylabel 'Residual'
title 'Scatter plot of the final residuals against the datum number'
