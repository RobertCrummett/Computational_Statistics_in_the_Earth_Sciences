%% Example 9.11
% Apply logistic regression the the iris flower data of Example 9.10 using
% the MATLAB function glmfit().

%% Importing Data
clear; close all; clf;
load fisheriris.mat

%% Data
x = [meas(51:150, :) ones(100, 1)];
y = [ones(1, 50) zeros(1, 50)]';

%% Logistic Regression
beta = glmfit(x(:, 1:4), y, "binomial", "link", "logit");

% Recalling that column one is the same as column five from Example 9.10,
% the results are nearly identical.