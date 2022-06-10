%% Example 4.14
% Compute 10000 random draws from the standardized Cauchy distribution and
% produce q-q plots of them against the standardized normal and Cauchy
% distributions. What differences do you see?

% The q-q plots are very useful for studying light-tailed data and
% detecting outliers, but is less useful for systematically heavy tailed
% data.

%% Random Draws
clear; close all; clf;
% MATLAB also provides qqplot() function. This here is a general formula to
% accomplish the same thing.
draws = 10000;
x = trnd(1,1,draws);
n = length(x);
u = ((1:n) - 0.5)/n;

%% Plotting q-q plots
% Normal quantiles
figure(1);
subplot(1,2,1); 
plot(norminv(u), sort(x), 'm+')
axis square
xlabel 'Normal Quantiles'; ylabel 'Sorted Data';

%Cauchy qunatiles
subplot(1,2,2); 
plot(tinv(u,1),sort(x),'m+')
axis square
xlabel 'Cauchy Quantiles';

% While the Gaussian distribution has exponential tails, the Cauchy
% distribution has algebraic tails.
% The sense of curvature in the q-q plots is opposite to that in the p-p
% plots.
% For heavy tailed data (Cauchy), the q-q plot has most of the data
% clustered in the middle and curvature controlled by outliers.