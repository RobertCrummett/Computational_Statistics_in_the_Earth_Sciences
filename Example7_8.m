%% Example 7.8
% The K-S test will be applied to a set of random draws from a standardized
% Gaussian distribution before and after it is contaminated with four
% outliers at the upper end of teh distribution.

%% Data
clear; close all; clf;
rng default;
data = normrnd(0, 1, 1000, 1);

%% K-S Test before Contamination
[h, p, ksstat, cv] = kstest(data);

%% Contaminate Data
figure(1);
qqplot(data)

data([1000, 999, 998, 997]) = [6.5 6.0 5.5 5];
figure(2);
qqplot(data)

% The dearture from normal is very apparent in figure 2.

%% K-S Test after Contamination
[h2, p2, ksstat2, cv2] = kstest(data);

% This demonstrates the insensitivity of the K-S test in the tails of the
% distribution. The p-value went up!