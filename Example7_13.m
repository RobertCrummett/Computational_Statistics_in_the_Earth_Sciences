%% Example 7.13
% Reanalyze the Gaussian simulation in Example 7.8 using the Jarque-Bera
% test without and with the four outliers.

% The Jarque-Bera test or skewness-kurtosis score is a test of the joint
% null hypothesis that the skewness and the excess kurtosis are zero
% against the alternative that they are not.

%% Data
clear; close all; clf;
rng default;
data = normrnd(0, 1, 1000, 1);

%% J-B Test without Contamination
[h, p, jbstat, cv] = jbtest(data);

%% Data Contamination
data([1000, 999, 998, 997]) = [6.5 6.0 5.5 5];

%% J-B Test with Contamination
[h2, p2, jbstat2, cv2] = jbtest(data);

% The test easily detected the presence of outliers.