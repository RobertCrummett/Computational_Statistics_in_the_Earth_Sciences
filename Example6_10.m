%% Example 6.10 - 11
% In 1879 and 1882, A. A. Michelson measured the speed of light. The
% measurements are presumed to be independant. Test the data t osee whether
% the sample means are different, and estimate the power of the test.

%% Importing Data
clear; close all; clf;
michelson1 = importdata("michelson1.dat");
michelson2 = importdata("michelson2.dat");

%% Plotting Gaussian q-q Plots

% Plotting By Hand
n1 = length(michelson1);
n2 = length(michelson2);
u1 = ((1:n1) - 0.5)/n1;
u2 = ((1:n2) - 0.5)/n2;

figure(1);
subplot(1,2,1);
plot(norminv(u1),sort(michelson1),'b+')
axis square
xlabel 'Normal Qunatiles'; ylabel 'Sorted Data';
subplot(1,2,2);
plot(norminv(u2),sort(michelson2),'b+')
axis square
xlabel 'Normal Qunatiles';

% MATLAB qqplot() Function
figure(2);
subplot(1,2,1);
qqplot(michelson1)
axis square
subplot(1,2,2);
qqplot(michelson2)
axis square

% The stair stepping in the data is characteristic of numerical roundoff
% that sometimes happens when data units are changed.

%% Hypothesis Testing
% Null hypothesis is that the means of the experiments are the same. The
% alternative hypothesis is that the means are not the same.

xbar1 = mean(michelson1);
xbar2 = mean(michelson2);
sig1 = std(michelson1, 1); % biased sample standard deviation
sig2 = std(michelson2, 1); % biased sample standard deviation

that = (xbar1 - xbar2)/sqrt(((n1*sig1^2 + n2*sig2^2)/(n1 + n2 - 2))*(1/n1 + 1/n2));
cv = tinv(0.975, n1 + n2 - 2);

% Since the test statistic is greater than the critical value, the null
% hypothesis at alpha = 0.05 can be rejected

pval = 2 * (1 - tcdf(that, n1 + n2 - 2));

% The low p-value strongly supports the alternative hypothesis.

%% Power Curve Plotted

sp = sqrt((n1*sig1^2 + n2*sig2^2)/(n1 + n2 - 2));
tcrit = tinv(0.975, n1 + n2 - 2);
mudiff = -100:0.1:100;
effect = mudiff/sp;
delta = sqrt(n1*n2/(n1 + n2))*effect; %Non-centrality parameter.
power = 1 - nctcdf(tcrit, n1 + n2 - 2, delta) + nctcdf(-tcrit, n1 + n2 - 2, delta);

figure(3);
plot(effect, power, 'k', 'LineWidth', 2)
xlabel 'Effect Size'; 
ylabel 'Power';
title 'Power Curve of the Michelson Speed of Light Data';

% The observed effect size is about 1.03, which indicates that there is a
% 99% chance of rejecting the null hypothesis when it is false. There is
% therefore high confidence in the test.

%% What if the variance's of the samples is different?
% Apply the unequal variance version of the t test to the Michelson speed
% of light data.

[h, p] = ttest2(michelson1, michelson2, 'Vartype', 'unequal');

% The conclusion remains unchanged. There is either no disparity between
% the variances of the data sets or the unequal variance t test does not do
% a good job of correction.