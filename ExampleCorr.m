%% Example: Correlation Coefficent
% The distribution of the bivariate correlation coefficent obtained by
% normalizing the off-diagonal elements of the unbiased sample covariance
% matrix by the square root of the product of its diagonal elements was
% obtained by Fisher (1915).

% Given a pair of random variables, with the condition of mutual
% independance, correlation coefficent is easily calculated. Once assuming
% that the joint distribution is Gaussian, the correlation coefficent pdf
% can be calculated as follows.

%% Plotting Correlation Coefficent PDF with changing Population Correlation Coefficent
% The pdf for the correlation coefficent with 10 degrees of freedom and a
% population correlation coefficent of rho
clear; close all; clf;

r = -1:0.01:1;
N = 10; % degrees of freedom
rho = [0.2 0.4 0.6 0.8]; % population correlation coefficents
pdfCorrCoef = zeros([size(r) size(rho)]);

for j = 1:length(rho)
    for i = 1:length(r)
        pdfCorrCoef(i,j) = (N-2)*gamma(N-1)*(1-rho(j)^2)^((N-1)/2)*(1-r(i)^2)^((N-4)/2)*...
            Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho(j)*r(i))/2)/(sqrt(2*pi)*gamma(N-1/2)*...
            (1-rho(j)*r(i))^(N-3/2));
    end
end

figure(1);
hold on
for i = 1:length(rho)
    plot(r,pdfCorrCoef(:,i))
end
hold off
xlabel 'Correlation Coefficent'; ylabel 'Probability Density';
legend('\rho = 0.2','\rho = 0.4','\rho = 0.6','\rho = 0.8','Location','nw')

%% Plotting Correlation Coefficent PDF with changing Degrees of Freedom
% The pdf for the correlation coefficent with 0 population value with
% changing degrees of freedom

r = -1:0.01:1;
N = [6 10 30 100]; % degrees of freedom
rho = 0; % population correlation coefficents
pdfCorrCoef = zeros(size(r));

for j = 1:length(N)
    for i = 1:length(r)
        pdfCorrCoef(i,j) = (N(j)-2)*gamma(N(j)-1)*(1-rho^2)^((N(j)-1)/2)*(1-r(i)^2)^((N(j)-4)/2)*...
            Hypergeometric2f1(1/2, 1/2, N(j)-1/2, (1+rho*r(i))/2)/(sqrt(2*pi)*gamma(N(j)-1/2)*...
            (1-rho*r(i))^(N(j)-3/2));
    end
end

figure(2);
hold on
for i = 1:length(N)
    plot(r,pdfCorrCoef(:,i))
end
hold off
xlabel 'Correlation Coefficent'; ylabel 'Probability Density';
legend('N = 6','N = 10','N = 30','N = 100','Location','nw')

%% Central Moments
% Needs work... change variables here to caculate different moments and
% variances

r = -1:0.01:1;
N = 10;
rho = 0.2;

fun = @(r, rho, N) r*(N-2)*gamma(N-1)*(1-rho^2)^((N-1)/2)*(1-r^2)^((N-4)/2)*...
        Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho*r)/2)/(sqrt(2*pi)*gamma(N-1/2)*...
        (1-rho*r)^(N-3/2));

%numeric integration for expected value
expectedVal = integral(@(r) fun(r, rho, N), -1, 1, 'ArrayValued', true);

%numeric integration for variance
fun1 = @(r, rho, N, xm) (r - xm)^2*(N-2)*gamma(N-1)*(1-rho^2)^((N-1)/2)*...
    (1-r^2)^((N-4)/2)*Hypergeometric2f1(1/2, 1/2, N-1/2, (1+rho*r)/2)/(sqrt(2*pi)*...
    gamma(N-1/2)*(1-rho*r)^(N-3/2));

xm = integral(@(r) fun(r, rho, N), -1, 1, 'ArrayValued', true); % the expected value
Var = integral(@(r) fun1(r, rho, N, xm), -1, 1, 'ArrayValued', true);