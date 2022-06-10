%% Example 4.12
% Compute 1000 random draws from the standardized Cauchy distribution, and
% produce stabilized p-p plots of them against the standardized normal and
% Cauchy distributions.
% What differences are there?

% Based on the probability integral transform.
% Deviations from a straight line are suggestive of the nature of the
% distribution. 

%% Random draws
clear; close all; clf;

picks = 1000;
x = trnd(1, [1 picks]);
u = ((1:picks) - 0.5)/picks; % uniform distribution

%% Arcsine transform applied to p-p plot
% Leads to constant variance at all points. Due to Michael (1983).
% This is not the case for the regular p-p plots. They have the highest
% variance near the mode and the lowest variance near the tails without
% this transform. So the transform is superior and always recommended.

r = 2/pi*asin(sqrt(u));
s1 = 2/pi*asin(sqrt(normcdf(sort(x)))); %Gaussian
s2 = 2/pi*asin(sqrt(tcdf(sort(x), 1))); %Cauchy


%% Plotting variance stabilized p-p plots
% Gaussian
figure(1);
subplot(1,2,1); axis square
plot(r, s1, 'k+ ')
ylabel 'Transformed Data'; xlabel 'Sine Quantiles';
title 'Standard Gaussian'
grid on

% Cauchy
subplot(1,2,2); axis square
plot(r, s2, 'k+ ')
xlabel 'Sine Quantiles'; 
title 'Standard Cauchy' 
grid on

% Stabilized p-p plot against the Gaussian is long tailed, as evidenced by
% the sharp downturn of the distribution tails. The Cauchy distribution is
% a straight line, as expected.
