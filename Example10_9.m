%% Example 10.9
% Returning to the temperature-humidity-evaporation data from Example 10.7,
% compute and assess the canonical correlation among the data.

%% Importing Data
clear; clf; close all;
the = importdata("the.dat");

%% Data
x = the(:, 1:9);
y = the(:, 10:end);

%% Canonical Correlation
[a, b, r, u, v, stats] = canoncorr(x, y);

disp(stats.pF)
disp(stats.pChisq)

% The canonical correlations are high, and the p-values are all small,
% rejecting the null hypothesis that the data correlations are zero.

%% Plot Canonical Variable Scores
figure(1);
plot(u(:,1), v(:,1), 'r+', 'LineWidth', 1.5)
hold on
plot(u(:,2), v(:,2), 'bo', 'LineWidth', 1.5)
hold off
axis equal square

xlabel 'U'
ylabel 'V'
legend(strcat("R^2 = ", string(r(1))), strcat("R^2 = ", string(r(2))),'Location','nw')
title 'Canonical Variable Scores for Temperature-Humidity-Wind Speed Data'