%% Example 9.6
% Revisit the star data with the Huber robust algorithm.

%% Importing Data
clear; close all; clf;
star = importdata('star.dat');

%% Data
y = star(:,2);
x = [star(:,1) ones(size(y))];

%% Robust M-Estimator
[tstruct, rstruct, pstruct] = Huber(y, x);

% This is an example of M-estimator weakness. It does not sense any of the
% leverage points evident in the hat matrix q-q plot. This is because
% M-estimators detect influential points only based on the size of the
% residuals.

%% Q-Q Plots for the Weighted Residuals and Hat Matrix
i = 1:length(y);
u = (i - 1/2)/length(y);

% Weighted Residual Q-Q Plot
figure(1);
plot(norminv(u), sort(rstruct.w(:,end).*rstruct.res(:,end)), '+k', 'LineWidth', 2)

xlabel 'Normal Quantiles'
ylabel 'Residual Order Statistic'
title 'Normal Q-Q Plot for the Weighted Residuals'

wx = repmat(rstruct.w(:,end), 1, size(x,2)).*x;
h = wx*inv(wx'*wx)*wx';
hat = diag(h)*size(x,1)/size(x,2);

% Weighted Hat Matrix Q-Q Plot
figure(2);
plot(betainv(u, size(x,2) , size(x,1) - size(x,2)), sort(hat), '+b', 'LineWidth', 2)
xlim([0 0.15])
ylim([0 5])
xlabel 'Beta (2, 45) Quantiles'
ylabel 'Hat Matrix Order Statistic'
title 'Weighted Hat Matrix Diagonal Q-Q Plot'

% These Q-Q plots can be compared and contrasted with the original plots
% from Example 9.1