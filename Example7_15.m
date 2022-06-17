%% Example 7.15
% There are two groups: a test group that was treated with a drug and a
% control group that was not. The latter is the frist 21 rows of the data.
% Test the null hypothesis that the median value of the differences between
% the test and control groups is zero.

%% Importing Data
clear; close all; clf;
remiss = importdata('remiss.dat');
x = remiss(1:21,2);
y = remiss(22:end,2);

%% Sign Test
[p, h, stats] = signtest(x, y);
xihat = stats.sign;

% The null hypothesis is strongly rejected. The difference is clearly illustrated
% in the figure.

%% Plotting Data
figure(1);
plot(remiss(1:21,1),x,'.b','MarkerSize',15)
hold on
plot(remiss(22:end,1),y,'rx','MarkerSize',5,'LineWidth',2)
hold off

xlabel 'Patient Number'
ylabel 'Remission Time'
xlim([0 22])
ylim([0 36])
legend('Control Group','Test Group','Location','best')