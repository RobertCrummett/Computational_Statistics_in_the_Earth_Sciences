%% Example: Distribution of the Interquartile Range
% Rather than compute the leading coefficent, it is easier to compute the
% integral and normalize the result to unity.

% Numerical appraoches to integration are required to solve this, except in
% the simplest situations (uniform distribution).

%% Plotting IQR PDF's
clear; close all; clf;

figure(1);
Nlist = [4 12 48 100];
w = 0:0.01:5;

hold on
for i = 1:length(Nlist)
    N = Nlist(i);
    fun = @(x, w, N) normcdf(x).^(N/4 - 1).*normpdf(x).*(normcdf(x + w) - normcdf(x))...
        .^(N/2 - 1).*normpdf(x + w).*(1 - normcdf(x + w)).^(N/4);
    f = zeros(size(w));
    s = integral2(@(x, w) fun(x, w, N), -Inf, Inf, 0, Inf); % in order to normalize to unity
    
    for j = 1:length(w)
        f(j) = integral(@(x) fun(x, w(j), N), -Inf, Inf)/s;
    end

    plot(0:0.01:5,f,'LineWidth',1.5)
end
hold off

ylabel 'IQR PDF';
title 'Sampling Distribution for the IQR of a standard Gaussian Random Variable'
legend('N = 4','N = 12','N = 48','N = 100','Location','ne')

%% Expected Value of Interquartile Range
% Standard Gaussian random variable
N = 4;

fun = @(x, w, N) normcdf(x).^(N/4 - 1).*normpdf(x).*(normcdf(x + w) - normcdf(x))...
    .^(N/2 - 1).*normpdf(x + w).*(1 - normcdf(x + w)).^(N/4);
fun1 = @(x, w, N) w.*normcdf(x).^(N/4 - 1).*normpdf(x).*(normcdf(x + w) - normcdf(x))...
    .^(N/2 - 1).*normpdf(x + w).*(1 - normcdf(x + w)).^(N/4);
e = integral2(@(x, w) fun1(x, w, N), -Inf, Inf, 0, Inf, "AbsTol",1e-12,"RelTol",1e-8);
e = e/integral2(@(x, w) fun(x, w, N), -Inf, Inf, 0, Inf);

% The expected value is slightly higher than the population value of
% 1.3161, and it is unclear whether this represents bias or numerical
% roundoff.