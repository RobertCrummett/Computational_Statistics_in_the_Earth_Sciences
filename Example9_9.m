%% Example 9.9
% Apply ridge regression and the lasso to the wheat data from Example 9.4.
% Estimate the number of predictor variables after removing the leverage
% points identified earlier.

%% Importing Data
clear; close all; clf;
wheat = importdata('wheat.dat');

%% Data
y = wheat(:,1);
x = wheat(:,2:7);

% Removing Leverage Points
y([17, 25, 34]') = [];
x([17, 25, 34]',:) = [];

%% Ridge Regression
lambda = [10:-0.1:1.1, 1:-0.01:0.11, 0.1:-0.001:0.011, 0.01:-0.0001:0.0011, 0.001:-0.00001:0.0001];

% Apply ridge regression
for i = 1:length(lambda)
    b = ridge(y, x, lambda(i), 1);
    yhat = y - x*b;
    mse(i) = sum((yhat - y).^2)/length(y);
end

b = ridge(y, x, lambda, 0);

% Plot MSE against the regularization parameter
figure(1);
plot(log10(lambda), mse, 'k-', 'LineWidth', 1.75)
xlim([-4 1])
xlabel 'Log_{10}\lambda'
ylabel 'Mean Squared Error'
title 'MSE as a Function of Ridge Regression Shrinkage Parameter \lambda'

% Plot Coefficents against the regularization parameter
[~, i] = find(lambda == 0.3);
figure(2);
plot(log10(lambda), b(2:end,:), 'LineWidth', 1.75)
xlabel 'Log_{10}\lambda'
ylabel 'Coefficient'
title 'Regression Parameters as a Function of Ridge Regression Shrinkage Parameter \lambda'
legend('First','Second','Third','Fourth','Fifth','Sixth','Location','best')

% The first, second, and fifth regression coefficents are all very small
% for lambda approximately equal to 0.01, yeilding estimates for the
% remaining coefficients that are in good agreement with the results of
% Example 9.4

%% Lasso Regression
[ b, fitinfo] = lasso(x, y, 'Lambda', lambda);

% Plot MSE against the regularization parameter
figure(3);
plot(log10(lambda), fitinfo.MSE, 'k-', 'LineWidth', 1.75)
xlabel 'Log_{10}\lambda'
ylabel 'Mean Squared Error'
title 'MSE as a Function of Lasso Shrinkage Parameter \lambda'

% Plot Coefficents against the regularization parameter
figure(4);
plot(log10(lambda), b, 'LineWidth', 1.75)
xlabel 'Log_{10}\lambda'
ylabel 'Coefficient'
title 'Regression Parameters as a Function of Lasso Shrinkage Parameter \lambda'
legend('First','Second','Third','Fourth','Fifth','Sixth','Location','nw')

% The behavior of lasso is different than ridge regression. Lasso produces
% abrubt changes in the coefficents and eliminating those that are not
% significant.

