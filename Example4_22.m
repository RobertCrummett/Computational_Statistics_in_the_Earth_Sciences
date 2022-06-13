%% Example 4.22
% Find the expected value of the rth order statistic for the uniform
% distribution on [0 1].

% Order statistics are obtained by sorting a group of random variables into
% ascending order. This sorting introduces dependance into the data. See H.
% A. David and David & Nagaraja (2003) for more.

%% Plotting Order Statisitics from Uniform Distirbution [0 1]
clear; close all; clf;
% F(x), the cdf, of the uniform distribution is jsut F(x) = x. Thus f(x),
% the pdf, is f(x) = 1.
% All of the formulas in this script are derivable from that.

% range and population
x = 0:0.01:1;
N = 10;

figure(1);
r = 1;
plot(x,betapdf(x, r, N - r + 1),'b')
hold on
r = 3;
plot(x,betapdf(x, r, N - r + 1),'r')
r = 7;
plot(x,betapdf(x, r, N - r + 1),'g')
r = 10;
plot(x,betapdf(x, r, N - r + 1),'m')
hold off

ylabel 'Probability Density';
title 'PDF of order statistics, population = 10, parent distribution is uniform [0 1]'
legend('Order statistic r = 1','Order statistic r = 3','Order statistic r = 7',...
    'Order statistic r = 10','Location','north')

%% Expected Value of Order Statistics from Uniform Distribution
r = [1 2 3 4 5 6 7 8 9 10];
muhat = r./(N + 1); % expected values