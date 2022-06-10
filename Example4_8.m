%% Example 4.8
% Obtain 100 and 1000 random draws from the standard exponential
% distribution and plot the empirical cdf and pdf for each

%% Draws from exponential distribution
x100 = exprnd(1,100,1);
x1000 = exprnd(1,1000,1);
[f1,x1] = ecdf(x100); %100 draws
[f2,x2] = ecdf(x1000); %1000 draws

%% Emirical cdf
figure(1);
plot(x1,f1,'b-',x2,f2,'r-')
title 'Empircal CDF of exponential dist for 100 and 1000 draws'
xlim([0 6])
ylabel 'Cumulative Probability'
legend('N = 100','N = 1000','Location','best')
% much smoother for the higher number of draws
% empirical cdf is an unbiased estimator of the population cdf

%% Empirical pdf
%using ecdfhist() 
figure(2);
subplot(1,2,1); ecdfhist(f1,x1,20);
xlim([0 6])
subplot(1,2,2); ecdfhist(f2,x2,20);
xlim([0 6])
% smoother for larger sample size

%using histogram()
figure(3);
subplot(1,2,1); histogram(x1,'NumBins',20,'Normalization','pdf');
xlim([0 6])
subplot(1,2,2); histogram(x2,'NumBins',20,'Normalization','pdf');
xlim([0 6])
% smoother for larger sample size