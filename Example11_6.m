%% Example 11.6
% Evaluate the svd of the MAR data.

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

%% Singular Value Decomposition
g = Close(geomean(mar));
[u, s, v] = svd(Clr(Perturbation(mar, repmat(1./g, length(mar), 1))));

disp(diag(s).^2./sum(diag(s).^2));

% The first two eigenvalues explain ~92% of the variance.
% The last singular value is zero, reflecting the fact that the Clr results
% in a signular covariance matrix.

% The first 10 right singular vectors are:

disp(v(:, 1:10))