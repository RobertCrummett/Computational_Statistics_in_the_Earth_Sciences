%% Example: Student t's Distributions
% The t distribution applies when neither mean nor variance of the random
% variables is known  a priori, and hence, must be estimated from the data.

% The difference between Gaussian and t distributions is especially
% important for small sample sizes. The are indistinguishable when the
% degrees of freedom exceeds 30.

%% Plotting t distributions against Normal distribution
% nu is the degrees of freedom
clear; close all; clf;
x = -4:0.1:4;

% nu = 1
figure(1);
plot(x,tpdf(x,1),'b')
hold on
% nu = 2
plot(x,tpdf(x,2),'r')
% nu = 5
plot(x,tpdf(x,5),'g')
% nu = 10
plot(x,tpdf(x,10),'m')
% normal distribution
pd = makedist('Normal');
plot(x,pdf(pd,x),'k')
hold off

ylabel 'Probability Density';
xlim([-4 4])
legend('\nu= 1','\nu = 2','\nu = 5','\nu = 10','N(0,1)', 'Location','ne')
title("Student's t distributions with varying degrees of freedom vs Normal distribution")

%% Plotting Noncentral t Distributions
% lambda is noncentrality parameter
x = -10:0.1:10;

% lambda  = 1, nu = 5
figure(2);
plot(x, nctpdf(x,5,1),'b')
% lambda = 3, nu = 5
hold on
plot(x,nctpdf(x,5,3),'r')
% lambda = 5, nu = 5
plot(x,nctpdf(x,5,5),'g')
% central t distribution, nu = 5
plot(x,tpdf(x,5),'k')
hold off

ylabel 'Probability Density'
xlim([-10 10]); ylim([0 0.4]);
legend('\lambda = 1','\lambda = 3','\lambda = 5','Central t distribution','Location','nw')
title 'Noncentral t distributions with \nu = 5 compared to Central t distribution with \nu = 5'