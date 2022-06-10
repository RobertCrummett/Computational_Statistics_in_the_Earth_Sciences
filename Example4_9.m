%% Example 4.9
% Obtain 10000 random draws, half from N(0,1) and half from N(3,1), compute
% the kernel density estimator for the data, and then vary the bandwidth
% and compare with the actual distribution
clear; close all; clf;

%% Random Draws
x = normrnd(0,1,1,5000);
x = [x normrnd(3,1,1,5000)];

%% KDE
[f,xi] = ksdensity(x);

%% Plotting KDE 
figure(1);
plot(xi,f,'b-','LineWidth',2)
title 'Kernel Density Estimates for 5000 N(0,1) and 5000 N(3,1)'

%% Varying Bandwidth and Plotting
hold on
[f,xi] = ksdensity(x,'bandwidth',0.1);
plot(xi,f,'r-','LineWidth',2)
[f,xi] = ksdensity(x,'bandwidth',1);
plot(xi,f,'g-','LineWidth',2)
plot(xi, 0.5*(normpdf(xi,0,1) + normpdf(xi,3,1)),'k-','LineWidth',2)
hold off
legend('KDE','KDE bw = 0.1','KDE bw = 1','N(0,1) & N(3,1)','Location','nw')

% too large a bandwidth fails to show bimodality, while too small a
% bandwidth makes the pdf appear erratic
% the default bandwidth is quite good