%% Example: Chi Squared Distribution
% Plot Chi squared distribuion with nu = 1, 2, 3, and 5

% Chi squared distribution is the distribution of the sum of squares of
% random samples from the standard normal distribution. The chi squared
% distribution is equal to the gamma distribution with alpha = nu/2, and
% beta = 1/2.

% When the means of the normal distributions used to compute the chi
% squared distribution are not zero, the distribution becomes the non
% central chi squared distribution

%% Plotting central chi squared distirbutions
% nu is called the degrees of freedom
clear; close all; clf;
% nu = 1
figure(1);
x = 0:0.1:10; % assume that x is normally distributed
plot(x,chi2pdf(x,1),'b','LineWidth',2)
ylabel 'Chi Squared PDF'
xlim([0 10]); ylim([0 1])

% nu = 2
hold on
plot(x,chi2pdf(x,2),'r','LineWidth',2)

% nu = 3
plot(x,chi2pdf(x,3),'g','LineWidth',2)

% nu = 5
plot(x,chi2pdf(x,5),'m','LineWidth',2)
hold off
legend('\nu = 1','\nu = 2','\nu = 3','\nu = 5','Location','ne')

%% Plotting non-central chi squared distributions
% lambda is the non centrality parameter

% lambda = 1, nu = 3
figure(2);
x = 0:0.1:20; % assume that x is normally distributed
plot(x,ncx2pdf(x,3,1),'b','LineWidth',2)
ylabel ' Noncentral Chi Squared PDF'
xlim([0 20]); ylim([0 0.2])

% lambda = 3, nu = 3
hold on
plot(x,ncx2pdf(x,3,3),'r','LineWidth',2)

% lambda = 5, nu = 3
plot(x,ncx2pdf(x,3,5),'g','LineWidth',2)

% lambda = 10, nu = 3
plot(x,ncx2pdf(x,3,10),'m','LineWidth',2)
hold off
legend('\lambda = 1','\lambda = 3','\lambda = 5','\lambda = 10','Location','ne')
