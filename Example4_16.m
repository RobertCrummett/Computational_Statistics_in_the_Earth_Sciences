%% Example 4.16
% Suppose the distribution of a ratio of a standard Rayleigh and standard
% Gaussian process is of interest. Estimate the random variable with a
% KDE.

%% Random draws
n = 10000;
x = raylrnd(1, n, 1); % standard Rayleigh set
y = normrnd(0, 1, n, 1); % standard Gaussian set

z = x./y; % ratio of interest

%% Characterizing the PDF
ksdensity(z,'npoints',n) %10000 kernels when n = 10000
xlim([-75 75]) %truncated abscissa
ylabel ' Simulated PDF';
title 'Simulated Ratio of Rayleigh and Gaussian variables'

% The resulting pdf is very sharply peaked and heavy tailed. Even with
% 10000 kernel centers, the middle is not adequately characterized.

% This example illustrates how the ratio of random variables will often
% result in algebraic tails and infinite variance.
% This is because the inverted probability density function often possesses
% algebraic tails. For example, Gaussian, Student's t, and exponential. 
% See Lehmann & Schaffer (1988) for more information.