%% Example 6.15
% The F test applies when the distribution of the test statistic follows
% the F distribution under the null hypothesis. The F test is used to
% compare two variances from two independant samples for which the
% underlying distributions are Gaussian.
% Thus, it serves as a test for homogeneity of variance.
% But since it is sensitive to departures from normality, it should be
% avoided for this purpose in favor of Barlett's M test.

%% Importing Data
clear; close all; clf;
data1 = importdata('michelson1.dat');
data2 = importdata('michelson2.dat');

%% Hypothesis Testing
% Computing maximum and minimum standard deviations
s1 = std(data1, 1);
s2 = std(data2, 1);
smax = max([s1 s2]);
smin = min([s1 s2]);

if s1 == smax
    nmax = length(data1);
    nmin = length(data2);
else
    nmax = length(data2);
    nmin = length(data1);
end

% Test statistic
Fhat = smax^2 / smin^2;

% Critical Value
cv = finv(0.975, nmax - 1, nmin - 1);

% Since the test statistic is greater than the critical value, the null
% hypothesis can be rejected.

pval = 2*(1 - fcdf(Fhat, nmax - 1, nmin - 1)); % two-sided test
pval1 = 1 - fcdf(Fhat, nmax-1, nmin-1); % upper-tail test

% The small p-value is strong evidence for the alternative hypothesis.

%% Plotting Power Curve of the F Test
% Since the data fail the null hypothesis of equal variance, the t test
% from prior examples is in doubt.

% Two-sided test
xlo = finv(0.025, nmax-1, nmin-1);
xhi = finv(0.975, nmax-1, nmin-1);
xhi1 = finv(0.95, nmax-1,nmin-1);

varrat = 1:0.01:3;
delta = nmax*(varrat - 1);
power = 1 - ncfcdf(xhi, nmax-1, nmin-1, delta) + ncfcdf(xlo, nmax-1, nmin-1, delta);
power1 = 1 - ncfcdf(xhi1, nmax-1, nmin-1, delta);

figure(1);
plot(varrat, power,"Color",'r')
hold on
plot(varrat, power1, "Color",'c')
xlabel 'Variance Ratio'; 
ylabel 'Power'; 
title 'Power of the F Test for Variance Homogeneity';
legend('Two-Sided Test','Upper-Tail Test','Location','nw')

% For the observed variance ratio of 1.92, the power of the F test is only
% about 55%, so there is not a whole lot of certainty about the outcome of
% the test.