%% Example 11.2 & 11.3
% For the MAR data, evaluate the geometric mean.

%% Importing Data
clear; clf; close all;
mar = readmatrix("MAR.xls");

%% Data Cleaning
% Keep only the major oxides
mar = mar(:, 3:13);

% Keep only data will all of the major oxides
mar = [mar(1:4, :)', mar(22, :)', mar(47:88, :)']';

% Remove redundant outliers at 37-38 and 39-40
% Remove redundant data at 24-26 and 30-36
mar = [mar(1:25, :)', mar(27:35, :)', mar(41:47, :)']';

mar = Close(mar, 100);

%% Closed Geometric Mean
% This is the Aitchison geometry's equivalent of the classical center of
% the data, the mean.

disp(Close(geomean(mar), 100))

%% Example 11.3
% For the MAR data, compute and analyze the normalized variation matrix and
% total variance.

[T, totvar] = NVariation(mar);