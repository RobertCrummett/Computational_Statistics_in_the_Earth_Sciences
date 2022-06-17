%% Example 7.22
% For the Michelson speed of light data, test the null hypothesis that the
% data have the same dispersion.

% The Ansari-Bradley Test is a two-sample nonparametric test of dispersion
% for continuous populations having the same median and shape.

%% Importing Data
clear; close all; clf;
michelson1 = importdata('michelson1.dat');
michelson2 = importdata('michelson2.dat');

%% A-B Test
[h, p] = ansaribradley(michelson1, michelson2);

% The test rejects the null hypothesis strongly. But the medians of the tests are
% different - to fit with the assumtions of the test, they need to be the
% same.

%% Normalizing Data by Subtracting Median
michelson1 = importdata('michelson1.dat') - median(michelson1);
michelson2 = importdata('michelson2.dat') - median(michelson2);

%% A-B Test on Normalized Data
[hmed, pmed] = ansaribradley(michelson1, michelson2);

% The test accepted the null hypothesis strongly! This should be contrasted
% to the F test and Barlett's M test that respectively rejected the null
% hypothesis and barely accepted it.