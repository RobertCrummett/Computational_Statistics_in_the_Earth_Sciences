%% Example 4.6
% Obtain 1000 realizations of 1000 random draws from the uniform
% distribution on the interval [-1 1]. Plot empiracal pdf of 1, 10, 100,
% & 1000 point averages.
% Does the central limit theorem (CLT) apply?
% Repeat for the Cauchy distribution. Does CLT apply? Why?

% The CLT states that as a large iid sample drawn from any distribution
% with a population mean mu and a finite variance sigma squared will
% yeild a sample mean x bar that is approximately normal with mean mu
% and variance sigma squared/N, where N is sample size.

%% Uniform Draws and Plots
clear; close all; clf;
x = unifrnd(-1, 1, 1000, 1000); %draws

figure(1);
hold on
subplot(2,2,1); histogram(x(1,:))
title 'No average'
subplot(2,2,2); histogram(mean(x(1:10,:)))
title '10 Point Average'
subplot(2,2,3); histogram(mean(x(1:100,:)))
title '100 Point Average'
subplot(2,2,4); histogram(mean(x(1:1000,:)))
title '1000 Point Average'
hold off
% Looks more Gaussian as the number of variables being averaged rises.
% The variance also decreases.

%% Cauchy Distribution
%(t distrivution with one degree of freedom)

x =  trnd(1, 100, 100); %Cauchy distribution is t distribution with 1 dof

figure(2);
hold on
subplot(1,3,1); histogram(x(1,:));
title 'No average'
subplot(1,3,2); histogram(mean(x(1:10,:)));
title '10 Point Average'
subplot(1,3,3); histogram(mean(x)); 
title '100 Point Average'
hold off
% Peak is not increasingly concentrated as N rises. 
% Thus, the CLT does not apply. 
% This makes sense, because for Cauchy distribution mu and variance are undefined.